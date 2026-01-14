/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/* This program tests if PnetCDF can report correct file formats */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <unistd.h> /* getopt() */

#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

static
int test_io(const char *out_path, /* ignored */
            const char *in_path,
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    char filename[512];
    int err, nerrs=0, fmt, ncid;

    sprintf(filename,"%s/test_cdf.nc%d",in_path, format);
    // printf("%s at %d: input filename %s\n",basename(__FILE__),__LINE__,filename);

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid);
    if (format == NC_FORMAT_NETCDF4 && PNETCDF_DRIVER_NETCDF4 == 0) {
        EXP_ERR(NC_ENOTBUILT)
        return 0;
    }
    else {
        CHECK_ERR
        if (err != NC_NOERR) {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            printf("Error in %s at %d: rank %d failed to open file %s\n",
                   basename(__FILE__), __LINE__, rank, filename);
            return 1;
        }
    }

    /* test NULL argument */
    err = ncmpi_inq_format(ncid, NULL); CHECK_ERR

    err = ncmpi_inq_format(ncid, &fmt); CHECK_ERR

    if (fmt != format) {
        printf("Error in %s at %d: expect CDF-%d format for file %s but got %d\n",
               __FILE__,__LINE__,format,filename,fmt);
        nerrs++;
    }
    err = ncmpi_close(ncid); CHECK_ERR

    /* test NULL argument */
    err = ncmpi_inq_file_format(filename, NULL); CHECK_ERR

    err = ncmpi_inq_file_format(filename, &fmt); CHECK_ERR

    if (fmt != format) {
        printf("Error in %s at %d: expect CDF-%d format for file %s but got %d\n",
               __FILE__,__LINE__,format,filename,fmt);
        nerrs++;
    }

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(nc_formats) / sizeof(int);
    opt.formats  = nc_formats;
    opt.ina      = 0; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 0; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "inquiring file formats", opt, test_io);

    MPI_Finalize();

    return err;
}

