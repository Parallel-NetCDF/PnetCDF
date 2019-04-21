/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 *
 *  Check if arguments start, count, stride, and imap are properly ignored
 *  when get/put a scalar variable.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>


static int
tst_fmt(char *filename, int cmode)
{
    int err, nerrs=0, ncid, varid, buf;
    MPI_Offset start[1], count[1], stride[1], imap[1];

    /* Test NetCDF-4 first as ncvalidator checks only classic formats */
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER | cmode,
                       MPI_INFO_NULL, &ncid); CHECK_ERR

    /* define a scalar variable of integer type */
    err = ncmpi_def_var(ncid, "scalar_var", NC_INT, 0, NULL, &varid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    buf = 1;
    start[0] = 1;
    count[0] = 2;
    stride[0] = 2;
    imap[0] = 2;

    /* put */
    err = ncmpi_put_var1_int_all(ncid, varid, NULL,  &buf); CHECK_ERR
    err = ncmpi_put_var1_int_all(ncid, varid, start, &buf); CHECK_ERR

    err = ncmpi_put_vara_int_all(ncid, varid, start, count, &buf); CHECK_ERR
    err = ncmpi_put_vara_int_all(ncid, varid, NULL, count, &buf); CHECK_ERR
    err = ncmpi_put_vara_int_all(ncid, varid, start, NULL, &buf); CHECK_ERR
    err = ncmpi_put_vara_int_all(ncid, varid, NULL, NULL, &buf); CHECK_ERR

    err = ncmpi_put_vars_int_all(ncid, varid, start, count, stride, &buf); CHECK_ERR
    err = ncmpi_put_vars_int_all(ncid, varid, NULL, count, stride, &buf); CHECK_ERR
    err = ncmpi_put_vars_int_all(ncid, varid, start, NULL, stride, &buf); CHECK_ERR
    err = ncmpi_put_vars_int_all(ncid, varid, start, count, NULL, &buf); CHECK_ERR
    err = ncmpi_put_vars_int_all(ncid, varid, NULL, NULL, NULL, &buf); CHECK_ERR

    err = ncmpi_put_varm_int_all(ncid, varid, start, count, stride, imap, &buf); CHECK_ERR
    err = ncmpi_put_varm_int_all(ncid, varid, NULL, NULL, NULL, NULL, &buf); CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid); CHECK_ERR

    err = ncmpi_inq_varid(ncid, "scalar_var", &varid); CHECK_ERR

    /* get */
    err = ncmpi_get_var1_int_all(ncid, varid, NULL,  &buf); CHECK_ERR
    err = ncmpi_get_var1_int_all(ncid, varid, start, &buf); CHECK_ERR

    err = ncmpi_get_vara_int_all(ncid, varid, start, count, &buf); CHECK_ERR
    err = ncmpi_get_vara_int_all(ncid, varid, NULL, count, &buf); CHECK_ERR
    err = ncmpi_get_vara_int_all(ncid, varid, start, NULL, &buf); CHECK_ERR
    err = ncmpi_get_vara_int_all(ncid, varid, NULL, NULL, &buf); CHECK_ERR

    err = ncmpi_get_vars_int_all(ncid, varid, start, count, stride, &buf); CHECK_ERR
    err = ncmpi_get_vars_int_all(ncid, varid, NULL, count, stride, &buf); CHECK_ERR
    err = ncmpi_get_vars_int_all(ncid, varid, start, NULL, stride, &buf); CHECK_ERR
    err = ncmpi_get_vars_int_all(ncid, varid, start, count, NULL, &buf); CHECK_ERR
    err = ncmpi_get_vars_int_all(ncid, varid, NULL, NULL, NULL, &buf); CHECK_ERR

    err = ncmpi_get_varm_int_all(ncid, varid, start, count, stride, imap, &buf); CHECK_ERR
    err = ncmpi_get_varm_int_all(ncid, varid, NULL, NULL, NULL, NULL, &buf); CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR
    return nerrs;
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    char filename[256], *hint_value;
    int err, nerrs=0, rank, nprocs, bb_enabled=0;

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

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for get/put scalar variables ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

#ifdef DEBUG
    if (nprocs > 1 && rank == 0)
        printf("Warning: %s is designed to run on 1 process\n", argv[0]);
#endif

    /* check whether burst buffering is enabled */
    if (inq_env_hint("nc_burst_buf", &hint_value)) {
        if (strcasecmp(hint_value, "enable") == 0) bb_enabled = 1;
        free(hint_value);
    }

    nerrs += tst_fmt(filename, 0);
    nerrs += tst_fmt(filename, NC_64BIT_OFFSET);
    if (!bb_enabled) {
#ifdef ENABLE_NETCDF4
        nerrs += tst_fmt(filename, NC_NETCDF4);
        nerrs += tst_fmt(filename, NC_NETCDF4 | NC_CLASSIC_MODEL);
#endif
    }
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

