/*********************************************************************
 *
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example shows how to use the chunking and compression features of
 * PnetCDF to write a 3D record variable of integer data type in parallel. It
 * first defines a netCDF variable of size time * global_ny * global_nx where
 *    global_ny == NY and
 *    global_nx == (NX * number of MPI processes).
 * The data partitioning pattern is a column-wise partitioning across all
 * processes. Each process writes a subarray of size ny * nx per record.
 *
 *    To compile:
 *        mpicc -O2 chunk_compress.c -o chunk_compress -lpnetcdf
 *
 * Example commands for MPI run and outputs from running ncmpidump on the
 * output netCDF file produced by this example program:
 *
 *    % mpiexec -n 4 ./chunk_compress /pvfs2/wkliao/testfile.nc
 *
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *            time = UNLIMITED ; // (0 currently)  <-- Not used anymore
 *            Y = 10 ;
 *            X = 16 ;
 *            _datablock_dim_0 = 65593 ;
 *            _datablock_dim_1 = 57 ;
 *    variables:
 *            int var ;
 *                    var:_ndim = 3 ;               <-- true ndims
 *                    var:_dimids = 0, 1, 2 ;       <-- dim IDs
 *                    var:_datatype = 4 ;           <-- true data type
 *                    var:_varkind = 1 ;            <-- compressed or not
 *                    var:_chunkdim = 1, 10, 4 ;    <-- chunk sizes
 *                    var:_filter = 2 ;             <-- filter ID
 *                    var:_metaoffset = 4LL ;       <-- offset to metadata block
 *            byte _datablock_0(_datablock_dim_0) ;
 *                    _datablock_0:_varkind = 2 ;
 *            byte _datablock_1(_datablock_dim_1) ;
 *                    _datablock_1:_varkind = 2 ;
 *
 *    // global attributes:
 *                    :_comressed = 1 ;  <-- file contains compressed variables or not
 *                    :_nwrite = 2 ;
 *                    :_recsize = 2LL ;  <-- true number of records
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <mpi.h>
#include <pnetcdf.h>

#define NY 10
#define NX 4

static int verbose;

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] [-k format] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [-k format] file format: 1 for CDF-1, 2 for CDF-2, 5 for CDF-5\n"
    "       [filename] output netCDF file name\n";
    fprintf(stderr, help, argv0);
}

/*----< pnetcdf_check_mem_usage() >------------------------------------------*/
/* check PnetCDF library internal memory usage */
static int
pnetcdf_check_mem_usage(MPI_Comm comm)
{
    int err, nerrs=0, rank;
    MPI_Offset malloc_size, sum_size;

    MPI_Comm_rank(comm, &rank);

    /* print info about PnetCDF internal malloc usage */
    err = ncmpi_inq_malloc_max_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && verbose)
            printf("maximum heap memory allocated by PnetCDF internally is %lld bytes\n",
                   sum_size);

        /* check if there is any PnetCDF internal malloc residue */
        err = ncmpi_inq_malloc_size(&malloc_size);
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }
    else if (err != NC_ENOTENABLED) {
        printf("Error at %s:%d: %s\n", __FILE__,__LINE__,ncmpi_strerror(err));
        nerrs++;
    }
    return nerrs;
}

/*----< pnetcdf_io() >-------------------------------------------------------*/
static int
pnetcdf_io(MPI_Comm comm, char *filename, int cmode)
{
    int i, j, rank, nprocs, err, nerrs=0;
    int ncid, varid, dimid[3], buf[NY][NX];
    MPI_Offset  global_ny, global_nx;
    MPI_Offset start[3], count[3];
    MPI_Info info;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_chunking", "enable");
    MPI_Info_set(info, "nc_chunk_default_filter", "zlib");

    /* create a new file for writing ----------------------------------------*/
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, info, &ncid); ERR
    MPI_Info_free(&info);

    /* the global array is NY * (NX * nprocs) */
    global_ny = NY;
    global_nx = NX * nprocs;

    for (i=0; i<NY; i++)
        for (j=0; j<NX; j++)
             buf[i][j] = rank;

    /* define dimensions x and y */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "Y",    global_ny,    &dimid[1]); ERR
    err = ncmpi_def_dim(ncid, "X",    global_nx,    &dimid[2]); ERR

    /* define a 2D variable of integer type */
    err = ncmpi_def_var(ncid, "var", NC_INT, 3, dimid, &varid); ERR

    /* set chunking (1st dimension should always be 1 for record variable) */
    int chunk_dim[3] = {1, NY, NX};
    err = ncmpi_var_set_chunk(ncid, varid, chunk_dim);; ERR

    /* use default filter */
    err = ncmpi_var_set_filter(ncid, varid, NC_FILTER_DEFLATE); ERR

    /* exit define mode */
    err = ncmpi_enddef(ncid); ERR

    /* set subarray start and count */
    start[0] = 0;
    start[1] = 0;
    start[2] = NX * rank;
    count[0] = 1;
    count[1] = NY;
    count[2] = NX;

    /* write to the 1st record */
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, &buf[0][0]); ERR

    /* write to the 2nd record */
    start[0] = 1;
    err = ncmpi_put_vara_int_all(ncid, varid, start, count, &buf[0][0]); ERR

    /* check the current record dimension size */
    MPI_Offset dim_len;
    err = ncmpi_inq_dimlen(ncid, 0, &dim_len); ERR
    if (rank == 0)
        printf("Time dimension length = %lld\n", dim_len);

    err = ncmpi_close(ncid); ERR

    return nerrs;
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    char filename[256];
    int i, rank, kind=0, cmode=0, nerrs=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 1;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hqk:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'k': kind = atoi(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (verbose && rank == 0) printf("%s: example of using put_vara APIs\n",__FILE__);

    switch (kind) {
        case(2): cmode = NC_64BIT_OFFSET; break;
        case(5): cmode = NC_64BIT_DATA;   break;
        default: cmode = 0;
    }

    nerrs += pnetcdf_io(MPI_COMM_WORLD, filename, cmode);

    nerrs += pnetcdf_check_mem_usage(MPI_COMM_WORLD);

    MPI_Finalize();
    return (nerrs > 0);
}

