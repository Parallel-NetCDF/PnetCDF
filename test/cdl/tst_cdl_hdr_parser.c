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

/*----< usage() >------------------------------------------------------------*/
static void usage (char *argv0) {
    char *help = "Usage: %s [OPTION] FILE\n\
       [-h] Print this help message\n\
       [-v] Verbose mode\n\
       [-o path] Output netCDF file path\n\
       FILE: Input CDL file path\n";
    fprintf (stderr, help, argv0);
}

int main(int argc, char **argv)
{
    extern int optind;
    char *outfile, *infile;
    int i, j, rank, err=0, nerrs=0, hid, verbose;

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

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for CDL header parser ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    err = cdl_hdr_open(infile, &hid);
    if (err != NC_NOERR) exit(1);

    if (verbose) printf("==================================================\n");

    /* create a new netcdf file */
    int ncid, cmode, format;

    err = cdl_hdr_inq_format(hid, &format); CHECK_ERR

    if (verbose) printf("CDF file format: CDF-%d\n", format);

    cmode = NC_CLOBBER;
    if (format == 2) cmode |= NC_64BIT_OFFSET;
    else if (format == 5) cmode |= NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, outfile, cmode, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    char *name;
    int ndims, *dimids, nvars, nattrs;
    void *value;
    nc_type xtype;
    MPI_Offset size, nelems;

    /* define dimensions */
    err = cdl_hdr_inq_ndims(hid, &ndims); CHECK_ERR
    if (verbose) printf("dim: ndims %d\n", ndims);

    for (i=0; i<ndims; i++) {
        int dimid;

        err = cdl_hdr_inq_dim(hid, i, &name, &size); CHECK_ERR
        if (verbose) printf("\t name %s size %lld\n",name, size);

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
                    printf("\t\tattr %s type %d nelems %lld (%s)\n",
                            name, xtype,nelems,(char*)value);
                else
                    printf("\t\tattr %s type %d nelems %lld\n",
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
                printf("\t name %s type %d nelems %lld (%s)\n",
                        name, xtype, nelems,(char*)value);
            else
                printf("\t name %s type %d nelems %lld\n",
                        name, xtype, nelems);
        }

        err = ncmpi_put_att(ncid, NC_GLOBAL, name, xtype, nelems, value); CHECK_ERR
    }
    err = ncmpi_close(ncid); CHECK_ERR

    err = cdl_hdr_close(hid); CHECK_ERR

    if (infile != NULL) free(infile);
    if (outfile != NULL) free(outfile);

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

