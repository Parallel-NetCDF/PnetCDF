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
    // printf("%s at %d: input filename %s\n",__FILE__,__LINE__,filename);

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid);
    if (format == NC_FORMAT_NETCDF4 && PNETCDF_DRIVER_NETCDF4 == 0) {
        EXP_ERR(NC_ENOTBUILT)
        return 0;
    }
    else
        CHECK_ERR

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

#if 0
void
tst_read_usage(char *argv0)
{
    char *base_name = basename(argv0);
    char *help =
    "Usage: %s [OPTIONS]...[filename]\n"
    "       [-h] Print help\n"
    "       [-q] quiet mode\n"
    "       [filename]: input netCDF folder name (default:.)\n";
    fprintf(stderr, help, base_name);
}

int tst_read(int         argc,
             char      **argv,
             char       *msg,
             loop_opts   opt,
             int       (*tst_body)(const char*, int, int, MPI_Info))
{
    extern int optind;
    extern char *optarg;
    char path[512];

    /* IDs for the netCDF file, dimensions, and variables. */
    int nprocs, rank, err, nerrs=0, quiet, coll_io;
    int i, a, d, r, m, b;

    MPI_Info info=MPI_INFO_NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    quiet = 0;

    while ((i = getopt(argc, argv, "hq")) != EOF)
        switch(i) {
            case 'q':
                quiet = 1;
                break;
            case 'h':
            default:  if (rank==0) tst_read_usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }

    if (argv[optind] == NULL) strcpy(path, ".");
    else                      snprintf(path, 256, "%s", argv[optind]);

    if (rank == 0) {
        char *cmd_str = (char *)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for %s", basename(argv[0]), msg);
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    MPI_Info_create(&info);

    for (i=0; i<opt.num_fmts; i++) {

        for (a=0; a<opt.num_ina; a++) {
        for (d=0; d<opt.num_drv; d++) {
        for (r=0; r<opt.num_ind; r++) {
        for (m=0; m<opt.num_chk; m++) {
        for (b=0; b<opt.num_bb;  b++) {

            if (a == 0) {
                MPI_Info_set(info, "nc_num_aggrs_per_node", "0");
            } else {
                MPI_Info_set(info, "nc_num_aggrs_per_node", "2");
            }

            if (d == 0) {
                MPI_Info_set(info, "nc_pncio", "enable");
            } else {
                MPI_Info_set(info, "nc_pncio", "disable");
            }

            if (r == 0) {
                MPI_Info_set(info, "romio_no_indep_rw", "false");
            } else {
                MPI_Info_set(info, "romio_no_indep_rw", "true");
            }

            if (m == 0) {
                MPI_Info_set(info, "nc_data_move_chunk_size", "1048576");
            } else {
                MPI_Info_set(info, "nc_data_move_chunk_size", "100");
            }

            if (b == 0) {
                MPI_Info_set(info, "nc_burst_buf", "disable");
            }
            else {
#ifdef ENABLE_BURST_BUFFER
                MPI_Info_set(info, "nc_burst_buf", "enable");
                MPI_Info_set(info, "nc_burst_buf_dirname", TESTOUTDIR);
                MPI_Info_set(info, "nc_burst_buf_overwrite", "enable");
#else
                continue;
#endif
            }

            for (coll_io=0; coll_io<opt.num_mod; coll_io++) {
                if (!quiet && rank == 0)
                    printf("a=%d d=%d r=%d m=%d b=%d coll_io=%d\n",
                            a,d,r,m,b,coll_io);

                /* NetCDF4 does not allow to extend number of record numbers in
                 * independent data mode. NC_ECANTEXTEND will be returned.
                 */
                if (opt.formats[i] == NC_FORMAT_NETCDF4 && coll_io == 0)
                    continue;

                nerrs = tst_body(path, opt.formats[i], coll_io, info);
                if (nerrs != NC_NOERR) goto err_out;
            }
        } /* loop b */
        } /* loop m */
        } /* loop r */
        } /* loop d */
        } /* loop a */
    }
    MPI_Info_free(&info);

    /* check if there is any malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has "OFFFMT" bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
    }

err_out:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs)
            printf(FAIL_STR, nerrs);
        else
            printf(PASS_STR);
    }

    MPI_Finalize();

    return (nerrs > 0);
}
#endif

int main(int argc, char **argv) {

#ifdef ENABLE_NETCDF4
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_NETCDF4, NC_FORMAT_64BIT_DATA};
#else
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};
#endif

    loop_opts opt;

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 0;
    opt.drv      = 1;
    opt.ind      = 1;
    opt.chk      = 0;
    opt.bb       = 0;
    opt.mod      = 1;
    opt.hdr_diff = 0;
    opt.var_diff = 0;

    return tst_main(argc, argv, "inquiring file formats", opt, test_io);
}

