/*
 *  Copyright (C) 2023, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests conflicted in-memory data types against external data
 * types and see if the correct error codes were returned.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o check_type check_type.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 check_type testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h> /* open(), lseek() */
#include <sys/stat.h>  /* open() */
#include <fcntl.h>     /* open() */
#include <unistd.h>    /* lseek() */
#include <string.h>
#include <strings.h>   /* strcasecmp() */
#include <libgen.h>    /* basename() */
#include <errno.h>    /* basename() */

#include <pnetcdf.h>

#include <testutils.h>

#define ATTR_EXP_ERR(expect_err,name) { \
    if (err != expect_err) { \
        printf("Error at line %d in %s: attr=%s err=%s\n", \
               __LINE__,__FILE__,name,ncmpi_strerrno(err)); \
        nerrs++; \
    }  \
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    char *attr_name[12] = {"attr_NC_NAT",
                           "attr_NC_BYTE",
                           "attr_NC_CHAR",
                           "attr_NC_SHORT",
                           "attr_NC_INT",
                           "attr_NC_FLOAT",
                           "attr_NC_DOUBLE",
                           "attr_NC_UBYTE",
                           "attr_NC_USHORT",
                           "attr_NC_UINT",
                           "attr_NC_INT64",
                           "attr_NC_UINT64"};
    char buf[1024];
    int i, j, err, nerrs=0, rank, ncid, max_type;
    MPI_Offset header_size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (format == NC_FORMAT_CLASSIC ||
        format == NC_FORMAT_64BIT_OFFSET ||
        format == NC_FORMAT_NETCDF4_CLASSIC)
        max_type = NC_DOUBLE;
    else
        max_type = NC_UINT64;

    for (i=0; i<1024; i++) buf[i]=0;

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    for (i=NC_BYTE; i<=max_type; i++) {
        /* create a new file (or truncate it to 0 length) */
        err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid); CHECK_ERR

        if (i == NC_CHAR) {
            for (j=0; j<3; j++) buf[j]='a'+j;
            err = ncmpi_put_att_text(ncid, NC_GLOBAL, attr_name[i], 3, buf); CHECK_ERR
        }
        else {
            err = ncmpi_put_att(ncid, NC_GLOBAL, attr_name[i], i, 3, buf); CHECK_ERR
        }
        ATTR_EXP_ERR(NC_NOERR, attr_name[i])

        err = ncmpi_close(ncid); CHECK_ERR

        /* reopen the file */
        err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_WRITE, info, &ncid); CHECK_ERR
        err = ncmpi_redef(ncid); CHECK_ERR

        err = ncmpi_del_att(ncid, NC_GLOBAL, attr_name[i]);
        ATTR_EXP_ERR(NC_NOERR, attr_name[i])

        /* call enddef to recalculate the header size */
        err = ncmpi_enddef(ncid); CHECK_ERR

        if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
            /* obtained updated header size */
            err = ncmpi_inq_header_size(ncid, &header_size); CHECK_ERR
        }

        err = ncmpi_close(ncid); CHECK_ERR

        /* This deletion of attribute will make the file size larger than
         * expected, i.e. larger than if no attribute is ever created. It is
         * not a fatal error. The file is still a valid NetCDF file. We should
         * expect running ncvalidator to print a warning message, such as
         * "file size (80) is larger than expected (48)!" See the fix in PR 99.
         */

        if (rank == 0) {
            off_t file_size;

            /* remove file type prefix substring */
            char *fname = remove_file_system_type_prefix(out_path);

            int fd = open(fname, O_RDONLY, 0666);

            if (fd == -1) {
                printf("Error: file open %s (%s)\n",fname,strerror(errno));
                return 1;
            }

            /* obtain file size */
            file_size = lseek(fd, 0, SEEK_END);

            if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC &&
                file_size != header_size)
                printf("Warning: expected file size "OFFFMT" but got %lld\n",
                       header_size, (long long)file_size);

            close(fd);
        }
    }

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
    opt.ind      = 0; /* test hint romio_no_indep_rw */
    opt.chk      = 1; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "delete attr", opt, test_io);

    MPI_Finalize();

    return err;
}
