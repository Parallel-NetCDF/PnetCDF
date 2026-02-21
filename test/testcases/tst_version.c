/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  Check whether PnetCDF version string returned from ncmpi_inq_libvers()
 *  matches the constant PNETCDF_VERSION defined in header file pnetcdf.h.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* strcpy(), strtok(), strcmp() */
#include <libgen.h>  /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

static
int test_io(const char *out_path, /* ignored */
            const char *in_path,  /* ignored */
            int         format,   /* ignored */
            int         coll_io,  /* ignored */
            MPI_Info    info)     /* ignored */
{
    char *str, *pnetcdf_version_str;
    int nerrs=0;

    str = (char*) malloc(strlen(ncmpi_inq_libvers())+1);
    strcpy(str, ncmpi_inq_libvers());
    pnetcdf_version_str = strtok(str, " ");

    if (pnetcdf_version_str == NULL) {
        printf("\nError: ncmpi_inq_libvers() returns ill form string %s\n",
               ncmpi_inq_libvers());
        nerrs++;
    }
    else if (strcmp(pnetcdf_version_str, PNETCDF_VERSION)) {
        printf("\nError: ncmpi_inq_libvers() returns %s does not match with PNETCDF_VERSION %s\n",
               pnetcdf_version_str, PNETCDF_VERSION);
        nerrs++;
    }
    free(str);

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {0};

    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 0; /* test intra-node aggregation */
    opt.drv      = 0; /* test PNCIO driver */
    opt.ind      = 0; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 0; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "PnetCDF library version", opt, test_io);

    MPI_Finalize();

    return err;
}
