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
int check_modes(char *filename)
{
    int rank, err, nerrs=0, ncid, cmode;
    char *path;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* delete the file and ignore error */
    /* remove the file system type prefix name if there is any.
     * For example, when filename = "lustre:/home/foo/testfile.nc", remove
     * "lustre:" to make path = "/home/foo/testfile.nc" in open() below
     */
    path = strchr(filename, ':');
    if (path == NULL) path = filename; /* no prefix */
    else              path++;

    if (rank == 0) unlink(path);
    MPI_Barrier(MPI_COMM_WORLD);

    /* create a new file and test various cmodes ----------------------------*/

    /* It is illegal to use both NC_64BIT_OFFSET and NC_64BIT_DATA together */
    cmode = NC_CLOBBER | NC_64BIT_OFFSET | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    EXP_ERR(NC_EINVAL_CMODE)

    /* The file should not be created */
    if (rank == 0) {
        if (access(path, F_OK) == 0) {
            printf("Error at %s:%d : file (%s) should not be created\n",
                   __FILE__,__LINE__, filename);
            nerrs++;
            /* delete the file and ignore error */
            unlink(path);
        }
        /* else : file does not exist */
    }
    MPI_Barrier(MPI_COMM_WORLD);

#ifdef ENABLE_NETCDF4
    /* It is illegal to use both NC_64BIT_OFFSET and NC_NETCDF4 together */
    cmode = NC_CLOBBER | NC_64BIT_OFFSET | NC_NETCDF4;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    EXP_ERR(NC_EINVAL_CMODE)

    /* The file should not be created */
    if (rank == 0) {
        if (access(path, F_OK) == 0) {
            printf("Error at %s:%d : file (%s) should not be created\n",
                   __FILE__,__LINE__, filename);
            nerrs++;
            /* delete the file and ignore error */
            unlink(path);
        }
        /* else : file does not exist */
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /* It is illegal to use both NC_64BIT_DATA and NC_NETCDF4 together */
    cmode = NC_CLOBBER | NC_64BIT_DATA | NC_NETCDF4;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    EXP_ERR(NC_EINVAL_CMODE)

    /* The file should not be created */
    if (rank == 0) {
        if (access(path, F_OK) == 0) {
            printf("Error at %s:%d : file (%s) should not be created\n",
                   __FILE__,__LINE__, filename);
            nerrs++;
            /* delete the file and ignore error */
            unlink(path);
        }
        /* else : file does not exist */
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    /* Collectively opening a non-existing file for read, expect error code
     * NC_ENOENT on all processes */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);

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
                       __LINE__,__FILE__, filename);
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
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_WRITE, MPI_INFO_NULL, &ncid);

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
                       __LINE__,__FILE__, filename);
                nerrs++;
                /* delete the file and ignore error */
                unlink(path);
            }
            /* else : file does not exist */
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid);
    CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

int main(int argc, char** argv)
{
    char *filename=NULL;
    int len, rank, err, nerrs=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) filename = strdup(argv[1]);
    else           filename = strdup("testfile.nc");
    len = strlen(filename) + 1;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(filename, len, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for file create/open modes ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* test under safe mode enabled */
    setenv("PNETCDF_SAFE_MODE", "1", 1);
    nerrs += check_modes(filename);

    /* test under safe mode disabled */
    setenv("PNETCDF_SAFE_MODE", "0", 1);
    nerrs += check_modes(filename);

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
    free(filename);

    MPI_Finalize();
    return (nerrs > 0);
}

