/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests if the MPI info object set inside the test program can be
 * merged with the ones set in the environment variable PNETCDF_HINTS. This
 * program was first developed to fix an internal bug in combine_env_hints()
 * where strings are allocated using strdup() but freed with NCI_Free.
 * See r3514.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_info tst_info.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 ./tst_info testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h> /* putenv() */
#include <string.h> /* strerror() */
#include <strings.h> /* strcasecmp() */
#include <errno.h>  /* errno */
#include <libgen.h> /* basename() */
#include <assert.h>
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define _RNDUP(x, unit) ((((x) + (unit) - 1) / (unit)) * (unit))

#define CHECK_HINT(hint) { \
    MPI_Info_get(info_used, hint, len, value, &flag); \
    if (!flag) { \
        printf("Error: hint \"%s\" is missing\n", hint); \
        nerrs++; \
    } \
}

static
int check_pnetcdf_hints(int ncid)
{
    char value[MPI_MAX_INFO_VAL+1];
    int err, nerrs=0, len=MPI_MAX_INFO_VAL, flag;
    MPI_Info info_used;

    err = ncmpi_inq_file_info(ncid, &info_used);
    CHECK_ERR;

    CHECK_HINT("nc_header_align_size");
    CHECK_HINT("nc_var_align_size");
    CHECK_HINT("nc_header_read_chunk_size");
    CHECK_HINT("nc_record_align_size");
    CHECK_HINT("pnetcdf_subfiling");
    CHECK_HINT("nc_in_place_swap");
    CHECK_HINT("nc_ibuf_size");

    MPI_Info_free(&info_used);

    return nerrs;
}

int main(int argc, char** argv) {
    char filename[256], value[MPI_MAX_INFO_VAL];
    int ncid1, ncid2, rank, err, nerrs=0, len, flag, varid;
    MPI_Offset header_size, header_extent, expect;
    MPI_Info info, info_used;

    MPI_Init(&argc, &argv);
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
        sprintf(cmd_str, "*** TESTING C   %s for merging env info ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* create a new file and keep it opened to make ncmpii_mem_root not NULL */
    err = ncmpi_create(MPI_COMM_WORLD, "dummy", NC_CLOBBER, MPI_INFO_NULL, &ncid1); CHECK_ERR

    /* retrieve MPI info object and check if all PnetCDF recognizable hints are
     * present */
    nerrs += check_pnetcdf_hints(ncid1);

    /* set environment variable PNETCDF_HINTS */
    err = setenv("PNETCDF_HINTS",
                 "romio_ds_write=disable;pnetcdf_subfiling=enable", 1);
    if (err != 0) {
        fprintf(stderr,"Error at line %d of %s calling setenv(): %s\n",
                __LINE__,__FILE__,strerror(errno));
        nerrs++;
    }

    /* create some PnetCDF-level I/O hints */
    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_header_align_size", "1");   /* size in bytes */
    MPI_Info_set(info, "nc_var_align_size",    "197"); /* size in bytes */

    /* create another new file using a non-NULL MPI info --------------------*/
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid2); CHECK_ERR

    MPI_Info_free(&info);

    err = ncmpi_def_var(ncid1, "var", NC_INT, 0, NULL, &varid); CHECK_ERR
    err = ncmpi_def_var(ncid2, "var", NC_INT, 0, NULL, &varid); CHECK_ERR

    /* set fill mode, so ncmpidiff can compare 2 output files without error */
    err = ncmpi_set_fill(ncid1, NC_FILL, NULL); CHECK_ERR
    err = ncmpi_set_fill(ncid2, NC_FILL, NULL); CHECK_ERR

    /* calling ncmpi_enddef() to write the file header */
    err = ncmpi_enddef(ncid2); CHECK_ERR

    /* NULL argument test */
    err = ncmpi_inq_header_size(ncid2, NULL); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid2, NULL); CHECK_ERR
    err = ncmpi_inq_file_info(ncid2, NULL); CHECK_ERR

    err = ncmpi_inq_header_size(ncid2, &header_size); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid2, &header_extent); CHECK_ERR

    /* get all the hints used by ncid2 */
    err = ncmpi_inq_file_info(ncid2, &info_used); CHECK_ERR

    err = ncmpi_close(ncid2); CHECK_ERR
    err = ncmpi_close(ncid1); CHECK_ERR

    /* delete file dummy */
    MPI_Barrier(MPI_COMM_WORLD);
    if (err == NC_NOERR && rank == 0) {
        err = ncmpi_delete("dummy", MPI_INFO_NULL); CHECK_ERR
    }

#ifdef VERBOSE
    if (rank == 0) {
        printf("\nheader_size = %lld\n",header_size);
        printf("header_extent=%lld\n\n",header_extent);
    }
#endif

    MPI_Info_get_valuelen(info_used, "nc_header_align_size", &len, &flag);
    if (flag) {
        MPI_Info_get(info_used, "nc_header_align_size", len+1, value, &flag);
        expect = _RNDUP(1, 4);
        if (expect != strtoll(value,NULL,10)) {
            printf("Error: nc_header_align_size expect %lld but got %lld\n",
                    expect, strtoll(value,NULL,10));
            nerrs++;
        }
    } else {
        printf("Error: hint \"nc_header_align_size\" is missing\n");
        nerrs++;
    }

    MPI_Info_get_valuelen(info_used, "nc_var_align_size", &len, &flag);
    if (flag) {
        MPI_Info_get(info_used, "nc_var_align_size", len+1, value, &flag);
        expect = _RNDUP(197, 4);
        if (expect != strtoll(value,NULL,10)) {
            printf("Error: nc_var_align_size expect %lld but got %lld\n",
                   expect, strtoll(value,NULL,10));
            nerrs++;
        }
    } else {
        printf("Error: hint \"nc_var_align_size\" is missing\n");
        nerrs++;
    }

    MPI_Info_get_valuelen(info_used, "romio_ds_write", &len, &flag);
    if (flag) {
        MPI_Info_get(info_used, "romio_ds_write", len+1, value, &flag);
        if (strcasecmp("disable", value)) {
            printf("Error: romio_ds_write expect \"disable\" but got \"%s\"\n",
                   value);
            nerrs++;
        }
    }

    MPI_Info_get_valuelen(info_used, "pnetcdf_subfiling", &len, &flag);
    if (flag) {
        MPI_Info_get(info_used, "pnetcdf_subfiling", len+1, value, &flag);
#ifdef ENABLE_SUBFILING
        if (strcasecmp("enable", value)) {
            printf("Error: pnetcdf_subfiling expect \"enable\" but got \"%s\"\n",
                   value);
            nerrs++;
        }
#else
        if (strcasecmp("disable", value)) {
            printf("Error: pnetcdf_subfiling expect \"disable\" but got \"%s\"\n",
                   value);
            nerrs++;
        }
#endif
    } else {
        printf("Error: hint \"pnetcdf_subfiling\" is missing\n");
        nerrs++;
    }
    MPI_Info_free(&info_used);

    /* set environment variable PNETCDF_HINTS */
    err = unsetenv("PNETCDF_HINTS");
    if (err != 0) {
        fprintf(stderr,"Error at line %d of %s calling unsetenv(): %s\n",
                __LINE__,__FILE__,strerror(errno));
        nerrs++;
    }

    /* re-open the file and get the MPI info object */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid1); CHECK_ERR

    /* retrieve MPI info object and check if all PnetCDF recognizable hints are
     * present */
    nerrs += check_pnetcdf_hints(ncid1);

    err = ncmpi_close(ncid1); CHECK_ERR

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

