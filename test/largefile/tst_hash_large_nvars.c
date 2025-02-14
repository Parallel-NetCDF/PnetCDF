/*********************************************************************
 *
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program is to test hash function performance using a large number of
 * variables e.g. > 100K
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

#define NVARS 400000

int main(int argc, char** argv)
{
    char filename[256];
    int i, rank, nprocs, err, nerrs=0, ncid, cmode, dimid, *varid, verbose=1;
    double timing[4], max_timing[4];
#ifdef PNC_MALLOC_TRACE
    MPI_Offset malloc_size[2], sum_size, max_size[2];
#endif
    MPI_Info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* get command-line arguments */
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
        sprintf(cmd_str, "*** TESTING C   %s for hashing of large nvars ",
                basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    if (verbose && rank == 0) printf("\nNVARS = %d\n", NVARS);

    varid = (int*) malloc(sizeof(int) * NVARS);

    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_hash_size_var", "2048");
    MPI_Info_set(info, "nc_hash_size_vattr", "2");

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    MPI_Info_free(&info);

#ifdef PNC_MALLOC_TRACE
    err = ncmpi_inq_malloc_size(&malloc_size[0]); CHECK_ERR
    err = ncmpi_inq_malloc_max_size(&malloc_size[1]); CHECK_ERR
    MPI_Reduce(&malloc_size, &max_size, 2, MPI_OFFSET, MPI_MAX, 0, MPI_COMM_WORLD);
    if (verbose && rank == 0) {
        printf("After ncmpi_create,  PnetCDF memory footprint high watermark %6.1f MB\n",
               (float)max_size[1]/1048576);
        printf("After ncmpi_create,  PnetCDF memory footprint                %6.1f MB\n",
               (float)max_size[0]/1048576);
    }
    fflush(stdout);
#endif

    err = ncmpi_def_dim(ncid, "dim", nprocs, &dimid);
    CHECK_ERR

    MPI_Barrier(MPI_COMM_WORLD);
    timing[0] = MPI_Wtime();
    for (i=0; i<NVARS; i++) {
        char name[64];
        long prefix, suffix;
        prefix = ((long)i*1747)%8642+100000;
        suffix = ((long)i*8313)%97531+100000;
        sprintf(name, "v%ld.x%ld", prefix, suffix);
        err = ncmpi_def_var(ncid, name, NC_INT, 1, &dimid, &varid[i]);
        CHECK_ERR
    }
    timing[0] = MPI_Wtime() - timing[0];

#ifdef PNC_MALLOC_TRACE
    err = ncmpi_inq_malloc_size(&malloc_size[0]); CHECK_ERR
    err = ncmpi_inq_malloc_max_size(&malloc_size[1]); CHECK_ERR
    MPI_Reduce(&malloc_size, &max_size, 2, MPI_OFFSET, MPI_MAX, 0, MPI_COMM_WORLD);
    if (verbose && rank == 0) {
        printf("After ncmpi_def_var, PnetCDF memory footprint high watermark %6.1f MB\n",
               (float)max_size[1]/1048576);
        printf("After ncmpi_def_var, PnetCDF memory footprint                %6.1f MB\n",
               (float)max_size[0]/1048576);
    }
    fflush(stdout);
#endif

    MPI_Barrier(MPI_COMM_WORLD);
    timing[1] = MPI_Wtime();
    for (i=0; i<NVARS; i++) {
        err = ncmpi_put_att(ncid, varid[i], "attr1", NC_INT, 1, &i);
        CHECK_ERR
        err = ncmpi_put_att(ncid, varid[i], "attr2", NC_INT, 1, &i);
        CHECK_ERR
    }
    timing[1] = MPI_Wtime() - timing[1];

#ifdef PNC_MALLOC_TRACE
    err = ncmpi_inq_malloc_size(&malloc_size[0]); CHECK_ERR
    err = ncmpi_inq_malloc_max_size(&malloc_size[1]); CHECK_ERR
    MPI_Reduce(&malloc_size, &max_size, 2, MPI_OFFSET, MPI_MAX, 0, MPI_COMM_WORLD);
    if (verbose && rank == 0) {
        printf("After ncmpi_put_att, PnetCDF memory footprint high watermark %6.1f MB\n",
               (float)max_size[1]/1048576);
        printf("After ncmpi_put_att, PnetCDF memory footprint                %6.1f MB\n",
               (float)max_size[0]/1048576);
    }
    fflush(stdout);
#endif

    MPI_Barrier(MPI_COMM_WORLD);
    timing[2] = MPI_Wtime();
    err = ncmpi_enddef(ncid);
    CHECK_ERR
    timing[2] = MPI_Wtime() - timing[2];

#ifdef PNC_MALLOC_TRACE
    err = ncmpi_inq_malloc_size(&malloc_size[0]); CHECK_ERR
    err = ncmpi_inq_malloc_max_size(&malloc_size[1]); CHECK_ERR
    MPI_Reduce(&malloc_size, &max_size, 2, MPI_OFFSET, MPI_MAX, 0, MPI_COMM_WORLD);
    if (verbose && rank == 0) {
        MPI_Offset header_size, header_extent;
        ncmpi_inq_header_extent(ncid, &header_extent);
        ncmpi_inq_header_size(ncid, &header_size);
        printf("After ncmpi_enddef,  PnetCDF memory footprint high watermark %6.1f MB\n",
               (float)max_size[1]/1048576);
        printf("After ncmpi_enddef,  PnetCDF memory footprint                %6.1f MB\n",
               (float)max_size[0]/1048576);
        printf("NetCDF file header size %lld extent %lld\n",header_size,header_extent);
    }
    fflush(stdout);
#endif

    /* Note the cost of ncmpi_close can be larger than expected when configured
     * with --enable-debug or --enable-profiling. This is because tracing
     * malloc builds a database and freeing a large number of memory pointers
     * involves searching and can be expensive.
     */
    MPI_Barrier(MPI_COMM_WORLD);
    timing[3] = MPI_Wtime();
    err = ncmpi_close(ncid);
    timing[3] = MPI_Wtime() - timing[3];
    CHECK_ERR

    free(varid);

#ifdef PNC_MALLOC_TRACE
    err = ncmpi_inq_malloc_size(&malloc_size[0]); CHECK_ERR
    err = ncmpi_inq_malloc_max_size(&malloc_size[1]); CHECK_ERR
    MPI_Reduce(&malloc_size, &max_size, 2, MPI_OFFSET, MPI_MAX, 0, MPI_COMM_WORLD);
    if (verbose && rank == 0) {
        printf("After ncmpi_close,   PnetCDF memory footprint high watermark %6.1f MB\n",
               (float)max_size[1]/1048576);
        printf("After ncmpi_close,   PnetCDF memory footprint                %4lld B\n",
               max_size[0]);
    }

    /* check if PnetCDF freed all internal malloc */
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size[0], &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0) {
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        }
    }
    if (malloc_size[0] > 0) ncmpi_inq_malloc_list();
#endif

    MPI_Reduce(&timing, &max_timing, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (verbose && rank == 0)
        printf("Time ncmpi_def_var = %.4f ncmpi_put_att = %.4f ncmpi_enddef = %.4f ncmpi_close = %.4f\n",
               max_timing[0],max_timing[1],max_timing[2],max_timing[3]);

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

