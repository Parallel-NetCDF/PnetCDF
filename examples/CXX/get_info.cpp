/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/*
 *    prints all MPI-IO hints used
 *    To compile:
 *        mpicxx -O2 get_info.c -o get_info -lpnetcdf
 *    To run:
 *        mpiexec -n 4 ./get_info [filename]
 *
 *    Example standard output:

    MPI File Info: nkeys = 18
    MPI File Info: [ 0] key =            cb_buffer_size, value = 16777216
    MPI File Info: [ 1] key =             romio_cb_read, value = automatic
    MPI File Info: [ 2] key =            romio_cb_write, value = automatic
    MPI File Info: [ 3] key =                  cb_nodes, value = 1
    MPI File Info: [ 4] key =         romio_no_indep_rw, value = false
    MPI File Info: [ 5] key =              romio_cb_pfr, value = disable
    MPI File Info: [ 6] key =         romio_cb_fr_types, value = aar
    MPI File Info: [ 7] key =     romio_cb_fr_alignment, value = 1
    MPI File Info: [ 8] key =     romio_cb_ds_threshold, value = 0
    MPI File Info: [ 9] key =         romio_cb_alltoall, value = automatic
    MPI File Info: [10] key =        ind_rd_buffer_size, value = 4194304
    MPI File Info: [11] key =        ind_wr_buffer_size, value = 524288
    MPI File Info: [12] key =             romio_ds_read, value = automatic
    MPI File Info: [13] key =            romio_ds_write, value = automatic
    MPI File Info: [14] key =            cb_config_list, value = *:1
    MPI File Info: [15] key =      nc_header_align_size, value = 512
    MPI File Info: [16] key =         nc_var_align_size, value = 512
    MPI File Info: [17] key = nc_header_read_chunk_size, value = 0
 */

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
using namespace std;

#include <string.h>
#include <pnetcdf>

using namespace PnetCDF;
using namespace PnetCDF::exceptions;

/*----< print_info() >------------------------------------------------------*/
static
void print_info(MPI_Info *info_used)
{
    int  i, nkeys;

    MPI_Info_get_nkeys(*info_used, &nkeys);
    printf("MPI File Info: nkeys = %d\n",nkeys);
    for (i=0; i<nkeys; i++) {
        char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];
        int  valuelen, flag;

        MPI_Info_get_nthkey(*info_used, i, key);
        MPI_Info_get_valuelen(*info_used, key, &valuelen, &flag);
        MPI_Info_get(*info_used, key, valuelen+1, value, &flag);
        printf("MPI File Info: [%2d] key = %25s, value = %s\n",i,key,value);
    }
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    char filename[128];
    int rank;
    MPI_Info info_used;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) strcpy(filename, argv[1]);
    else           strcpy(filename, "testfile.nc");

    try {
        /* create the file */
        NcmpiFile ncFile(MPI_COMM_WORLD, filename, NcmpiFile::replace,
                         NcmpiFile::data64bits);

        /* get all the default hints used */
        ncFile.Inq_file_info(&info_used);
    }
    catch(NcmpiException& e) {
       cout << e.what() << " error code=" << e.errorCode() << " Error!\n";
       return 1;
    }

    if (rank == 0) print_info(&info_used);
    MPI_Info_free(&info_used);

    MPI_Finalize();
    return 0;
}

