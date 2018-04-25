/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example sets two PnetCDF hints:
 *    nc_header_align_size and nc_var_align_size
 * and prints the hint values as well as the header size, header extent, and
 * two variables' starting file offsets.
 *
 * The compile and run commands are given below.
 *
 *    % mpicxx -O2 -o hints hints.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 ./hints /pvfs2/wkliao/testfile.nc
 *
 *    nc_header_align_size      set to = 1024
 *    nc_var_align_size         set to = 512
 *    nc_header_read_chunk_size set to = 256
 *    header size                      = 252
 *    header extent                    = 1024
 *    var_zy start file offset         = 1024
 *    var_yx start file offset         = 3072
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
using namespace std;

#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <pnetcdf>

using namespace PnetCDF;
using namespace PnetCDF::exceptions;

#define NZ 5
#define NY 5
#define NX 5

static void
usage(char *argv0)
{
    cerr <<
    "Usage: %s [-h] | [-q] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [filename] output netCDF file name\n"
    << argv0;
}

static
void print_hints(NcmpiFile &ncFile,
                 NcmpiVar  &var0,
                 NcmpiVar  &var1)
{
    char value[MPI_MAX_INFO_VAL];
    int len, flag;
    MPI_Offset header_size, header_extent, var_zy_start, var_yx_start;
    MPI_Offset h_align=-1, v_align=-1, h_chunk=-1;
    MPI_Info info_used;

    ncFile.Inq_header_size(&header_size);
    ncFile.Inq_header_extent(&header_extent);

    var0.Inq_file_offset(&var_zy_start);
    var1.Inq_file_offset(&var_yx_start);

    /* get all the hints used */
    ncFile.Inq_file_info(&info_used);

    MPI_Info_get_valuelen(info_used, (char*)"nc_header_align_size", &len, &flag);
    if (flag) {
        MPI_Info_get(info_used, (char*)"nc_header_align_size", len+1, value, &flag);
        h_align = strtoll(value,NULL,10);
    }
        MPI_Info_get_valuelen(info_used, (char*)"nc_var_align_size", &len, &flag);
    if (flag) {
        MPI_Info_get(info_used, (char*)"nc_var_align_size", len+1, value, &flag);
        v_align = strtoll(value,NULL,10);
    }
    MPI_Info_get_valuelen(info_used, (char*)"nc_header_read_chunk_size", &len, &flag);
    if (flag) {
        MPI_Info_get(info_used, (char*)"nc_header_read_chunk_size", len+1, value,&flag);
        h_chunk = strtoll(value,NULL,10);
    }
    MPI_Info_free(&info_used);

    if (h_align == -1)
        printf("nc_header_align_size      is NOT set\n");
    else
        printf("nc_header_align_size      set to = %lld\n", h_align);

    if (v_align == -1)
        printf("nc_var_align_size         is NOT set\n");
    else
        printf("nc_var_align_size         set to = %lld\n", v_align);
    if (h_chunk == -1)
        printf("nc_header_read_chunk_size is NOT set\n");
    else
        printf("nc_header_read_chunk_size set to = %lld\n", h_chunk);

    printf("header size                      = %lld\n", header_size);
    printf("header extent                    = %lld\n", header_extent);
    printf("var_zy start file offset         = %lld\n", var_zy_start);
    printf("var_yx start file offset         = %lld\n", var_yx_start);
}

int main(int argc, char** argv)
{
    extern int optind;
    char filename[256];
    int i, rank, nprocs, verbose=1, *buf_zy;
    float *buf_yx;
    vector<MPI_Offset> start(2), count(2);
    MPI_Info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hq")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

    try {
        MPI_Info_create(&info);
        MPI_Info_set(info, (char*)"nc_header_align_size",      (char*)"1024");
        MPI_Info_set(info, (char*)"nc_var_align_size",         (char*)"512");
        MPI_Info_set(info, (char*)"nc_header_read_chunk_size", (char*)"256");
        /* note that set the above values to 1 to disable the alignment */

        /* create a new file for writing -------------------------------------*/
        NcmpiFile ncFile(MPI_COMM_WORLD, filename, NcmpiFile::replace,
                         NcmpiFile::classic5, info);
        MPI_Info_free(&info);

        /* define 3 dimensions */
        vector<NcmpiDim> dimid(3);
        dimid[0] = ncFile.addDim("Z", NZ*nprocs);
        dimid[1] = ncFile.addDim("Y", NY*nprocs);
        dimid[2] = ncFile.addDim("X", NX*nprocs);

        /* define a variable of size (NZ * nprocs) * (NY * nprocs) */
        vector<NcmpiDim> dimid_zy(2);
        dimid_zy[0] = dimid[0];
        dimid_zy[1] = dimid[1];
        NcmpiVar var0 = ncFile.addVar("var_zy", ncmpiInt, dimid_zy);

        /* define a variable of size (NY * nprocs) * (NX * nprocs) */
        vector<NcmpiDim> dimid_yx(2);
        dimid_yx[0] = dimid[1];
        dimid_yx[1] = dimid[2];
        NcmpiVar var1 = ncFile.addVar("var_yx", ncmpiFloat, dimid_yx);

        /* var_zy is partitioned along Z dimension */
        buf_zy = (int*) malloc(NZ * (NY * nprocs) * sizeof(int));
        for (i=0; i<NZ*(NY*nprocs); i++) buf_zy[i] = i;

        start[0] = NZ * rank; start[1] = 0;
        count[0] = NZ;        count[1] = NY * nprocs;
        var0.putVar_all(start, count, buf_zy);

        /* var_yx is partitioned along X dimension */
        buf_yx = (float*) malloc((NY * nprocs) * NX * sizeof(float));
        for (i=0; i<(NY*nprocs)*NX; i++) buf_yx[i] = i;

        start[0] = 0;           start[1] = NX * rank;
        count[0] = NY * nprocs; count[1] = NX;
        var1.putVar_all(start, count, buf_yx);

        if (rank == 0 && verbose) print_hints(ncFile, var0, var1);

        free(buf_zy);
        free(buf_yx);
    }
    catch(NcmpiException& e) {
       cout << e.what() << " error code=" << e.errorCode() << " Error!\n";
       return 1;
    }

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    int err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Finalize();
    return 0;
}

