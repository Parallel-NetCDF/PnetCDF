/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This program tests if the correct error codes are returned given various
 * create/open modes.
 *
 * NC_EINVAL_CMODE should be returned when creating a file using
 * comde with both NC_64BIT_OFFSET & NC_64BIT_DATA flags set.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strdup() */
#include <libgen.h> /* basename() */
#include <unistd.h> /* unlink(), access() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define EXPECT_ERR2(err_no1, err_no2) \
    if (err != err_no1 && err != err_no2) { \
        nerrs++; \
        printf("Error at line %d in %s: expect error code %s or %s but got %s\n", \
               __LINE__,__FILE__,ncmpi_strerrno(err_no1),ncmpi_strerrno(err_no2),ncmpi_strerrno(err)); \
    }

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,  /* ignored */
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    int rank, err, nerrs=0, ncid, cmode;
    char *path;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* delete the file and ignore error */

    /* remove the file system type prefix name if there is any.  For example,
     * when out_path = "lustre:/home/foo/testfile.nc", remove "lustre:" to make
     * path pointing to "/home/foo/testfile.nc", so it can be used in POSIX
     * unlink() and access() below
     */
    path = remove_file_system_type_prefix(out_path);

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) unlink(path);

    MPI_Barrier(MPI_COMM_WORLD);

    /* create a new file and test various cmodes ----------------------------*/

    /* It is illegal to use both NC_64BIT_OFFSET and NC_64BIT_DATA together */
    cmode = NC_CLOBBER | NC_64BIT_OFFSET | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid);
    EXP_ERR(NC_EINVAL_CMODE)

    MPI_Barrier(MPI_COMM_WORLD);

    /* The file should not be created */
    if (rank == 0) {
        if (access(path, F_OK) == 0) {
            printf("Error at %s:%d : file (%s) should not be created\n",
                   __FILE__,__LINE__, out_path);
            nerrs++;
            /* delete the file and ignore error */
            unlink(path);
        }
        /* else : file does not exist */
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC) {
        /* It is illegal to use both NC_64BIT_OFFSET and NC_NETCDF4 together */
        if (format == NC_FORMAT_NETCDF4_CLASSIC)
            cmode = NC_CLOBBER | NC_64BIT_OFFSET | NC_NETCDF4 | NC_CLASSIC_MODEL;
        else
            cmode = NC_CLOBBER | NC_64BIT_OFFSET | NC_NETCDF4;

        err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid);
        EXP_ERR(NC_EINVAL_CMODE)

        MPI_Barrier(MPI_COMM_WORLD);

        /* The file should not be created */
        if (rank == 0) {
            if (access(path, F_OK) == 0) {
                printf("Error at %s:%d : file (%s) should not be created\n",
                       __FILE__,__LINE__, out_path);
                nerrs++;
                /* delete the file and ignore error */
                unlink(path);
            }
            /* else : file does not exist */
        }
        MPI_Barrier(MPI_COMM_WORLD);

        /* It is illegal to use both NC_64BIT_DATA and NC_NETCDF4 together */
        if (format == NC_FORMAT_NETCDF4_CLASSIC)
            cmode = NC_CLOBBER | NC_64BIT_DATA | NC_NETCDF4 | NC_CLASSIC_MODEL;
        else
            cmode = NC_CLOBBER | NC_64BIT_DATA | NC_NETCDF4;

        err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid);
        EXP_ERR(NC_EINVAL_CMODE)

        MPI_Barrier(MPI_COMM_WORLD);

        /* The file should not be created */
        if (rank == 0) {
            if (access(path, F_OK) == 0) {
                printf("Error at %s:%d : file (%s) should not be created\n",
                       __FILE__,__LINE__, out_path);
                nerrs++;
                /* delete the file and ignore error */
                unlink(path);
            }
            /* else : file does not exist */
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    /* Collectively opening a non-existing file for read, expect error code
     * NC_ENOENT on all processes */
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid);

    MPI_Barrier(MPI_COMM_WORLD);

    /* When using MVAPICH2 2.2, its Lustre driver adds O_CREAT to all open
     * calls. This is considered a bug in an MPI-IO implementation. Due to this
     * bug, the non-existing file will be created with zero-length and thus
     * PnetCDF spews NC_ENOTNC */
    if (err == NC_ENOTNC) {
        /* ignore the error and delete the file */
        if (rank == 0) unlink(path);
    }
    else {
        /* older version of OpenMPI and MPICH may return MPI_ERR_IO instead of
         * MPI_ERR_NO_SUCH_FILE */
        EXPECT_ERR2(NC_ENOENT, NC_EFILE)

        /* The file should not be created */
        if (rank == 0) {
            if (access(path, F_OK) == 0) {
                printf("Error at line %d in %s: file (%s) should not be created\n",
                       __LINE__,__FILE__, out_path);
                nerrs++;
                /* delete the file and ignore error */
                unlink(path);
            }
            /* else : file does not exist */
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /* Collectively opening a non-existing file for write, expect error code
     * NC_ENOENT on all processes */
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_WRITE, info, &ncid);

    MPI_Barrier(MPI_COMM_WORLD);

    /* When using MVAPICH2 2.2, its Lustre driver adds O_CREAT to all open
     * calls. This is considered a bug in an MPI-IO implementation. Due to this
     * bug, the non-existing file will be created with zero-length and thus
     * PnetCDF spews NC_ENOTNC */
    if (err == NC_ENOTNC) {
        /* ignore the error and delete the file */
        if (rank == 0) unlink(path);
    }
    else {
        /* older version of OpenMPI and MPICH may return MPI_ERR_IO instead of
         * MPI_ERR_NO_SUCH_FILE */
        EXPECT_ERR2(NC_ENOENT, NC_EFILE)

        /* The file should not be created */
        if (rank == 0) {
            if (access(path, F_OK) == 0) {
                printf("Error at line %d in %s: file (%s) should not be created\n",
                       __LINE__,__FILE__, out_path);
                nerrs++;
                /* delete the file and ignore error */
                unlink(path);
            }
            /* else : file does not exist */
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

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
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 0; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "file create/open modes", opt, test_io);

    MPI_Finalize();

    return err;
}
