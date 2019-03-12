/*********************************************************************
 *
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program shows how to use a single vard API call to write or read two
 * consecutive variables. This example uses two MPI datatype constructors,
 * MPI_Type_create_subarray and MPI_Type_create_hindexed, to create the same
 * subarray view for two variables, but with different lower and upper bound
 * MPI type maps. The two datatypes are then concatenated into a single
 * filetype.
 *
 *    To compile:
 *        mpicc -O2 vard_mvars.c -o vard_mvars -lpnetcdf
 *
 * Example commands for MPI run and outputs from running ncmpidump on the
 * NetCDF file produced by this example program:
 *
 *  % mpiexec -n 4 ./vard_mvars /pvfs2/wkliao/testfile.nc
 *
 * The expected results from the output file contents are:
 *
 *  % ncmpidump testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-1
 *    dimensions:
 *           REC_DIM = UNLIMITED ; // (2 currently)
 *           Y = 2 ;
 *           X = 20 ;
 *    variables:
 *           int fix_var0(Y, X) ;
 *           int fix_var1(Y, X) ;
 *           int rec_var2(REC_DIM, Y, X) ;
 *           int rec_var3(REC_DIM, Y, X) ;
 *    data:
 *
 *    fix_var0 =
 *      0, 1, 2, 3, 4, 100, 101, 102, 103, 104, 200, 201, 202, 203, 204, 300, 301,
 *       302, 303, 304,
 *     10, 11, 12, 13, 14, 110, 111, 112, 113, 114, 210, 211, 212, 213, 214, 310,
 *       311, 312, 313, 314 ;
 *
 *    fix_var1 =
 *      1000, 1001, 1002, 1003, 1004, 1100, 1101, 1102, 1103, 1104, 1200, 1201,
 *        1202, 1203, 1204, 1300, 1301, 1302, 1303, 1304,
 *      1010, 1011, 1012, 1013, 1014, 1110, 1111, 1112, 1113, 1114, 1210, 1211,
 *        1212, 1213, 1214, 1310, 1311, 1312, 1313, 1314 ;
 *
 *    rec_var2 =
 *      _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 *      _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 *      2000, 2001, 2002, 2003, 2004, 2100, 2101, 2102, 2103, 2104, 2200, 2201,
 *        2202, 2203, 2204, 2300, 2301, 2302, 2303, 2304,
 *      2010, 2011, 2012, 2013, 2014, 2110, 2111, 2112, 2113, 2114, 2210, 2211,
 *        2212, 2213, 2214, 2310, 2311, 2312, 2313, 2314 ;
 *
 *    rec_var3 =
 *      _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 *      _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 *      3000, 3001, 3002, 3003, 3004, 3100, 3101, 3102, 3103, 3104, 3200, 3201,
 *        3202, 3203, 3204, 3300, 3301, 3302, 3303, 3304,
 *      3010, 3011, 3012, 3013, 3014, 3110, 3111, 3112, 3113, 3114, 3210, 3211,
 *        3212, 3213, 3214, 3310, 3311, 3312, 3313, 3314 ;
 * }
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h> /* getopt() */

#include <mpi.h>
#include <pnetcdf.h>

#define NY 2
#define NX 5

static int verbose;

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

#define CHECK_VALUE(buf,base) { \
    for (j=0; j<NY; j++) { \
        for (i=0; i<NX; i++) \
            if ((buf)[j*NX+i] != (base)+rank*100+j*10+i) { \
                printf("line %d: expecting buf[%d*NX+%d]=%d but got %d\n",\
                       __LINE__,j,i,(base)+rank*100+j*10+i,(buf)[j*NX+i]); \
                nerrs++; \
            } \
    } \
}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
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

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    extern int optind;
    char         filename[256];
    int          i, j, err, ncid, varid[4], dimids[3], nerrs=0;
    int          rank, nprocs, *buf[2];
    int          array_of_sizes[2], array_of_subsizes[2], array_of_starts[2];
    int          array_of_blocklengths[NY];
    MPI_Offset   recsize, start[2], count[2], offset[2];
    MPI_Aint     a0, a1, array_of_displacements[NY];
    MPI_Datatype buftype, vtype[2], filetype;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 1;

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

    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (verbose && rank == 0) printf("%s: example of using vard APIs to write/read two variables\n",__FILE__);

    buf[0] = (int*)malloc(NY * NX * sizeof(int));
    for (j=0; j<NY; j++) for (i=0; i<NX; i++)
        buf[0][j*NX+i] = rank*100 + j*10 + i;

    buf[1] = (int*)malloc(NY * NX * sizeof(int));
    for (j=0; j<NY; j++) for (i=0; i<NX; i++)
        buf[1][j*NX+i] = 1000 + rank*100 + j*10 + i;

    /* construct buftype: concatenate two separated buffers */
    array_of_blocklengths[0] = NY*NX;
    array_of_blocklengths[1] = NY*NX;
    array_of_displacements[0] = 0;
    MPI_Get_address(buf[0], &a0);
    MPI_Get_address(buf[1], &a1);
    array_of_displacements[1] = a1 - a0;
    vtype[0] = vtype[1] = MPI_INT;
    MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
                           vtype, &buftype);
    MPI_Type_commit(&buftype);

    /* create a new NetCDF file */
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL,
                       &ncid); ERR

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "REC_DIM", NC_UNLIMITED, &dimids[0]); ERR
    err = ncmpi_def_dim(ncid, "Y",       NY,           &dimids[1]); ERR
    err = ncmpi_def_dim(ncid, "X",       NX*nprocs,    &dimids[2]); ERR

    /* define 2D fixed-size variables */
    err = ncmpi_def_var(ncid, "fix_var0", NC_INT, 2, dimids+1, &varid[0]); ERR
    err = ncmpi_def_var(ncid, "fix_var1", NC_INT, 2, dimids+1, &varid[1]); ERR
    /* define 3D record variables */
    err = ncmpi_def_var(ncid, "rec_var2", NC_INT, 3, dimids, &varid[2]); ERR
    err = ncmpi_def_var(ncid, "rec_var3", NC_INT, 3, dimids, &varid[3]); ERR

    /* enable fill mode for the two record variables */
    err = ncmpi_def_var_fill(ncid, varid[2], 0, NULL); ERR
    err = ncmpi_def_var_fill(ncid, varid[3], 0, NULL); ERR

    err = ncmpi_enddef(ncid); ERR

    /* specify the access pattern, a subarray for each variable of size
     * NY * NX from each MPI process */
    start[0] = 0;  start[1] = NX*rank;
    count[0] = NY; count[1] = NX;

    /* create the first datatype using MPI subarray constructor */
    array_of_sizes[0]    = 2;
    array_of_sizes[1]    = NX*nprocs;
    array_of_subsizes[0] = count[0];
    array_of_subsizes[1] = count[1];
    array_of_starts[0]   = start[0];
    array_of_starts[1]   = start[1];
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes,
                             array_of_starts, MPI_ORDER_C, MPI_INT, &vtype[0]);
    MPI_Type_commit(&vtype[0]);

    /* create the second datatype using MPI hindexed constructor */
    for (i=0; i<NY; i++) /* there are NY blocks, each of size NX */
        array_of_blocklengths[i] = NX;

    array_of_displacements[0] = NX*rank*sizeof(int);
    for (i=1; i<NY; i++)
        /* distance between any 2 consecutive blocks is NX*nprocs */
        array_of_displacements[i] = array_of_displacements[i-1]
                                  + NX*nprocs*sizeof(int);

    MPI_Type_create_hindexed(NY, array_of_blocklengths, array_of_displacements,
                             MPI_INT, &vtype[1]);
    MPI_Type_commit(&vtype[1]);

    /* concatenate two datatypes into filetype */
    array_of_blocklengths[0] = 1;
    array_of_blocklengths[1] = 1;
    array_of_displacements[0] = 0;

    /* calculate distance between 2 fixed-size variables */
    err = ncmpi_inq_varoffset(ncid, varid[0], &offset[0]); ERR
    err = ncmpi_inq_varoffset(ncid, varid[1], &offset[1]); ERR
    array_of_displacements[1] = offset[1] - offset[0];

    MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
                           vtype, &filetype);
    MPI_Type_commit(&filetype);

    /* write 2 consecutive fixed-size variables in a single vard call */
    err = ncmpi_put_vard_all(ncid, varid[0], filetype, buf[0], 1, buftype); ERR

    /* check if the contents of write buffers are altered */
    CHECK_VALUE(buf[0], 0)
    CHECK_VALUE(buf[1], 1000)

    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[0][j*NX+i] = -1;
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[1][j*NX+i] = -1;

    /* read back fixed-size variables */
    err = ncmpi_get_vard_all(ncid, varid[0], filetype, buf[0], 1, buftype); ERR

    /* check the contents of read buffers */
    CHECK_VALUE(buf[0], 0)
    CHECK_VALUE(buf[1], 1000)

    MPI_Type_free(&filetype);

    /* obtain record size (sum of individual record of all variables) */
    err = ncmpi_inq_recsize(ncid, &recsize); ERR

    /* calculate distance between 2 record variables */
    err = ncmpi_inq_varoffset(ncid, varid[2], &offset[0]); ERR
    err = ncmpi_inq_varoffset(ncid, varid[3], &offset[1]); ERR
    array_of_displacements[1] = offset[1] - offset[0];

    /* fill the 1st record of two variables with fill value */
    err = ncmpi_fill_var_rec(ncid, varid[2], 0); ERR
    err = ncmpi_fill_var_rec(ncid, varid[3], 0); ERR

    /* skip the 1st record, and write to the 2nd record of two consecutive
     * variables */
    array_of_displacements[0] += recsize;
    array_of_displacements[1] += recsize;

    MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
                           vtype, &filetype);
    MPI_Type_commit(&filetype);

    for (j=0; j<NY; j++) for (i=0; i<NX; i++)
        buf[0][j*NX+i] = 2000 + rank*100 + j*10 + i;

    for (j=0; j<NY; j++) for (i=0; i<NX; i++)
        buf[1][j*NX+i] = 3000 + rank*100 + j*10 + i;

    /* write 2 consecutive record variables in a single call to vard API */
    err = ncmpi_put_vard_all(ncid, varid[2], filetype, buf[0], 1, buftype); ERR

    /* check if the contents of write buffers are altered */
    CHECK_VALUE(buf[0], 2000)
    CHECK_VALUE(buf[1], 3000)

    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[0][j*NX+i] = -1;
    for (j=0; j<NY; j++) for (i=0; i<NX; i++) buf[1][j*NX+i] = -1;

    /* read back record variables */
    err = ncmpi_get_vard_all(ncid, varid[2], filetype, buf[0], 1, buftype); ERR

    /* check the contents of read buffers */
    CHECK_VALUE(buf[0], 2000)
    CHECK_VALUE(buf[1], 3000)

    MPI_Type_free(&vtype[0]);
    MPI_Type_free(&vtype[1]);
    MPI_Type_free(&filetype);
    MPI_Type_free(&buftype);
    free(buf[0]); free(buf[1]);

    err = ncmpi_close(ncid); ERR

    nerrs += pnetcdf_check_mem_usage(MPI_COMM_WORLD);

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    MPI_Finalize();
    return (nerrs > 0);
}
