/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program reads the header of a text file in CDL format and creates a new
 * netCDF file with the same header. This is to test the CDL header parser APIs.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_cdl_hdr_parser tst_cdl_hdr_parser.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 tst_cdl_hdr_parser -i cdl_header.txt -o testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> /* getopt() */
#include <libgen.h> /* basename() */

#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    char *name;
    int i, j, rank, err=0, nerrs=0, hid, ncid, verbose;
    int ndims, *dimids, nvars, nattrs;
    void *value;
    nc_type xtype;
    MPI_Offset size, nelems;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 0;

    err = cdl_hdr_open(in_path, &hid);
    if (err != NC_NOERR) exit(1);

    if (verbose) printf("==================================================\n");

    /* create a new netcdf file */
    err = cdl_hdr_inq_format(hid, &format); CHECK_ERR

    if (verbose) printf("CDF file format: CDF-%d\n", format);

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* define dimensions */
    err = cdl_hdr_inq_ndims(hid, &ndims); CHECK_ERR
    if (verbose) printf("dim: ndims %d\n", ndims);

    for (i=0; i<ndims; i++) {
        int dimid;

        err = cdl_hdr_inq_dim(hid, i, &name, &size); CHECK_ERR
        if (verbose) printf("\t name %s size "OFFFMT"\n",name, size);

        err = ncmpi_def_dim(ncid, name, size, &dimid); CHECK_ERR
    }

    /* define variables */
    err = cdl_hdr_inq_nvars(hid, &nvars); CHECK_ERR
    if (verbose) printf("var: nvars %d\n", nvars);

    for (i=0; i<nvars; i++) {
        int varid;

        err = cdl_hdr_inq_var(hid, i, &name, &xtype, &ndims, &dimids); CHECK_ERR

        err = ncmpi_def_var(ncid, name, xtype, ndims, dimids, &varid); CHECK_ERR

        /* define local attributes */
        err = cdl_hdr_inq_nattrs(hid, i, &nattrs); CHECK_ERR

        if (verbose) {
            printf("\t name %s type %d ndims %d nattr %d\n",
                          name, xtype, ndims, nattrs);
            for (j=0; j<ndims; j++)
                printf("\t\tdimid %d\n",dimids[j]);
        }

        for (j=0; j<nattrs; j++) {
            err = cdl_hdr_inq_attr(hid, i, j, &name, &xtype, &nelems, &value); CHECK_ERR
            if (verbose) {
                if (xtype == NC_CHAR)
                    printf("\t\tattr %s type %d nelems "OFFFMT" (%s)\n",
                            name, xtype,nelems,(char*)value);
                else
                    printf("\t\tattr %s type %d nelems "OFFFMT"\n",
                           name, xtype, nelems);
            }

            err = ncmpi_put_att(ncid, varid, name, xtype, nelems, value); CHECK_ERR
        }
    }

    /* define global attributes */
    err = cdl_hdr_inq_nattrs(hid, NC_GLOBAL, &nattrs); CHECK_ERR
    if (verbose) printf("global attrs: nattrs %d\n", nattrs);

    for (i=0; i<nattrs; i++) {
        err = cdl_hdr_inq_attr(hid, NC_GLOBAL, i, &name, &xtype, &nelems, &value);
        if (verbose) {
            if (xtype == NC_CHAR)
                printf("\t name %s type %d nelems "OFFFMT" (%s)\n",
                        name, xtype, nelems,(char*)value);
            else
                printf("\t name %s type %d nelems "OFFFMT"\n",
                        name, xtype, nelems);
        }

        err = ncmpi_put_att(ncid, NC_GLOBAL, name, xtype, nelems, value); CHECK_ERR
    }
    err = ncmpi_close(ncid); CHECK_ERR

    err = cdl_hdr_close(hid); CHECK_ERR

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(nc_formats) / sizeof(int);
    opt.formats  = nc_formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 1; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "CDL header parser", opt, test_io);

    MPI_Finalize();

    return err;
}

