/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This example program reads the header of a text file in CDL format and
 * creates a new netCDF file with the same metadata.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o create_from_cdl create_from_cdl.c -lpnetcdf
 *
 *    % mpiexec -n 4 create_from_cdl -i cdl_header.txt -o testfile.nc
 *
 * File header of the newly created file, testfile.nc, should be the same as
 * the input CDL file, except for the number of records (i.e. dimension
 * NC_UNLIMITED). This fact can be verified by running 'ncmpidump/ncdump'
 * utility with command-line option '-h'.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> /* getopt() */
#include <libgen.h> /* basename() */

#include <mpi.h>
#include <pnetcdf.h>

static int verbose;

#define CHECK_ERR { \
    if (err != NC_NOERR) { \
        printf("Error at %s:%d : %s\n", __FILE__,__LINE__, \
               ncmpi_strerror(err)); \
        nerrs++; \
    } \
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

/*----< usage() >------------------------------------------------------------*/
static void usage (char *argv0) {
    char *help = "Usage: %s [OPTION] FILE\n\
       [-h] Print this help message\n\
       [-v] Verbose mode\n\
       [-o path] Output netCDF file path\n\
       FILE: Input CDL file path (required)\n";
    fprintf (stderr, help, argv0);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv)
{
    extern int optind;
    char *outfile, *infile, *name;
    int i, j, rank, err=0, nerrs=0, hid;
    int ncid, cmode, format, ndims, nvars, nattrs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 1;
    infile = NULL;
    outfile = NULL;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hqo:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'o': outfile = strdup(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (outfile == NULL) outfile = strdup("testfile.nc");

    if (argv[optind] == NULL) {
        if (rank == 0) usage(argv[0]);
        MPI_Finalize();
        return 1;
    }
    infile = strdup(argv[optind]);

    /* read and parse the input CDL header file */
    err = cdl_hdr_open(infile, &hid);
    if (err != NC_NOERR) exit(1);

    /* obtain file format version */
    err = cdl_hdr_inq_format(hid, &format);
    CHECK_ERR

    /* create a new netcdf file */
    cmode = NC_CLOBBER;
    if (format == 2) cmode |= NC_64BIT_OFFSET;
    else if (format == 5) cmode |= NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, outfile, cmode, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    /* retrieve the number of dimensions defined in the CDL file */
    err = cdl_hdr_inq_ndims(hid, &ndims);
    CHECK_ERR
    if (verbose) printf("Number of dimensions : %d\n", ndims);

    for (i=0; i<ndims; i++) {
        int dimid;
        MPI_Offset dimlen;

        /* retrieve metadata of dimension i */
        err = cdl_hdr_inq_dim(hid, i, &name, &dimlen);
        CHECK_ERR

        /* define a new dimension in the new file */
        err = ncmpi_def_dim(ncid, name, dimlen, &dimid);
        CHECK_ERR
    }

    /* retrieve number of variables defined in the CDL file */
    err = cdl_hdr_inq_nvars(hid, &nvars);
    CHECK_ERR
    if (verbose) printf("Number of variables : %d\n", nvars);

    for (i=0; i<nvars; i++) {
        int varid, *dimids;
        nc_type xtype;

        /* retrieve metadata of variable i defined in the CDL file */
        err = cdl_hdr_inq_var(hid, i, &name, &xtype, &ndims, &dimids);
        CHECK_ERR

        /* define a new variable in the new file */
        err = ncmpi_def_var(ncid, name, xtype, ndims, dimids, &varid);
        CHECK_ERR

        /* retrieve metadata of attribute j associated with variable i */
        err = cdl_hdr_inq_nattrs(hid, i, &nattrs);
        CHECK_ERR

        for (j=0; j<nattrs; j++) {
            void *value;
            nc_type xtype;
            MPI_Offset nelems;

            /* retrieve metadata of attribute j associated with variable i */
            err = cdl_hdr_inq_attr(hid, i, j, &name, &xtype, &nelems, &value);
            CHECK_ERR

            /* define a new attribute and associate it to variable, varid */
            err = ncmpi_put_att(ncid, varid, name, xtype, nelems, value);
            CHECK_ERR
        }
    }

    /* retrieve the number of global attributes */
    err = cdl_hdr_inq_nattrs(hid, NC_GLOBAL, &nattrs);
    CHECK_ERR
    if (verbose) printf("Number of global attributes: %d\n", nattrs);

    for (i=0; i<nattrs; i++) {
        void *value;
        nc_type xtype;
        MPI_Offset nelems;

        /* retrieve metadata of global attribute i */
        err = cdl_hdr_inq_attr(hid, NC_GLOBAL, i, &name, &xtype, &nelems,
                               &value);
        CHECK_ERR

        /* define a global attribute in the new file */
        err = ncmpi_put_att(ncid, NC_GLOBAL, name, xtype, nelems, value);
        CHECK_ERR
    }

    /* close the NetCDF file */
    err = ncmpi_close(ncid);
    CHECK_ERR

    /* close the CDL file */
    err = cdl_hdr_close(hid);
    CHECK_ERR

    if (infile != NULL) free(infile);
    if (outfile != NULL) free(outfile);

    nerrs += pnetcdf_check_mem_usage(MPI_COMM_WORLD);

    MPI_Finalize();
    return (nerrs > 0);
}

