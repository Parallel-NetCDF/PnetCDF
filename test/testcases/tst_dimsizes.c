/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program (borrowed from netCDF library) tests defining and inquiring the
 * maximum allowable dimension size for CDF-1, 2, and 5 formats.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_dimsizes tst_dimsizes.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 tst_dimsizes testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

#define DIMMAXCLASSIC NC_MAX_INT
#define DIMMAX64OFFSET NC_MAX_INT
#define DIMMAX64DATA NC_MAX_INT64

/*
 * NetCDF file format specification:
 *   netcdf_file  = header data
 *   header       = magic numrecs dim_list gatt_list var_list
 *   dim_list     = ABSENT | NC_DIMENSION nelems [dim ...]
 *   dim          = name dim_length
 *   dim_length   = NON_NEG
 *   NON_NEG      = <non-negative INT>
 *   INT          = <32-bit signed integer, Bigendian, two's complement>
 *
 * Therefore, the max dimension size are:
 * NC_CLASSIC       Max dimension size is  NC_INT_MAX
 * NC_64BIT_OFFSET  Max dimension size is  NC_INT_MAX
 * NC_64BIT_DATA    Max dimension size is  NC_INT64_MAX
 * Note that for NC_64BIT_DATA, the max dimension size is different from netCDF
 * library. This is because PnetCDF uses MPI_Offset for dimension size and
 * MPI_Offset is a signed long long.
*/

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int err, nerrs=0, fmt, ncid, dimid;
    MPI_Offset dimsize;

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* Writing Max Dimension Size */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR
    dimsize = -1;
    err = ncmpi_def_dim(ncid, "testdim1", dimsize, &dimid); EXP_ERR(NC_EDIMSIZE)

    if (format == NC_FORMAT_CLASSIC) {
        dimsize = DIMMAXCLASSIC;
        err = ncmpi_def_dim(ncid, "testdim", dimsize, &dimid); CHECK_ERR
        dimsize = (MPI_Offset)DIMMAXCLASSIC+1;
        err = ncmpi_def_dim(ncid, "testdim1", dimsize, &dimid); EXP_ERR(NC_EDIMSIZE)
    }
    else if (format == NC_FORMAT_64BIT_OFFSET) {
        dimsize = DIMMAX64OFFSET;
        err = ncmpi_def_dim(ncid, "testdim", dimsize, &dimid); CHECK_ERR
        dimsize = (MPI_Offset)DIMMAX64OFFSET+1;
        err = ncmpi_def_dim(ncid, "testdim1", dimsize, &dimid); EXP_ERR(NC_EDIMSIZE)
    }
    else {
        dimsize = DIMMAX64DATA;
        err = ncmpi_def_dim(ncid, "testdim", dimsize, &dimid); CHECK_ERR
        dimsize = -1;
        err = ncmpi_def_dim(ncid, "testdim1", dimsize, &dimid); EXP_ERR(NC_EDIMSIZE)
    }

    err = ncmpi_close(ncid); CHECK_ERR

    /* Reading Max Dimension Size */
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOCLOBBER, info, &ncid); CHECK_ERR
    err = ncmpi_inq_format(ncid, &fmt); CHECK_ERR
    err = ncmpi_inq_dimid(ncid, "testdim", &dimid); CHECK_ERR
    err = ncmpi_inq_dimlen(ncid, dimid, &dimsize); CHECK_ERR
    if (format == NC_FORMAT_CLASSIC && dimsize != DIMMAXCLASSIC) {
        printf("Error at line %d in %s: expecting dimsize %d but got "OFFFMT"\n",
                __LINE__,__FILE__,DIMMAXCLASSIC,dimsize);
        nerrs++;
    }
    else if (format == NC_FORMAT_64BIT_OFFSET && dimsize != DIMMAX64OFFSET) {
        printf("Error at line %d in %s: expecting dimsize %d but got "OFFFMT"\n",
                __LINE__,__FILE__,DIMMAX64OFFSET,dimsize);
        nerrs++;
    }
    else if (format == NC_FORMAT_64BIT_DATA && dimsize != DIMMAX64DATA) {
        printf("Error at line %d in %s: expecting dimsize %lld but got "OFFFMT"\n",
                __LINE__,__FILE__,(long long)DIMMAX64DATA,dimsize);
        nerrs++;
    }

    err = ncmpi_close(ncid); CHECK_ERR

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
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "defining max dimension sizes", opt, test_io);

    MPI_Finalize();

    return err;
}
