/*
 *  Copyright (C) 2023, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/*
 * This program tests if ncmpi_create() can keep the file a symbolic link when
 * NC_CLOBBER is used and the symbolic link file is already exist.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>    /* strdup() */
#include <strings.h>   /* strcasecmp() */
#include <libgen.h>    /* basename() */
#include <errno.h>     /* errno */

#include <sys/types.h> /* lstat() */
#include <sys/stat.h>  /* lstat() */
#include <unistd.h>    /* lstat(), symlink(), unlink(), sync() */

#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define ERR_HANDLE(name, err) {                                      \
    if (err != 0) {                                                  \
        fprintf(stderr,"Error at file %s line %d: %s failed (%s)\n", \
                __FILE__,__LINE__, name, strerror(errno));           \
        nerrs++;                                                     \
    }                                                                \
}

int main(int argc, char **argv) {
    char *filename, *symlink_fname, *fname;
    int  err, nerrs=0, len, rank, ncid, verbose=0;
    struct stat statbuf;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) filename = strdup(argv[1]);
    else           filename = strdup("testfile.nc");
    len = (int)strlen(filename) + 1;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(filename, len, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for NC_CLOBBER on symlink file ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* remove file system type prefix substring */
    fname = remove_file_system_type_prefix(filename);

    symlink_fname = (char*) malloc(strlen(filename) + 10);

    /* create a regular file and a symbolic link to it */
    err = 0;
    if (rank == 0) {
        FILE *fp;

        /* when calling POSIX I/O, remove file type prefix from file name */
        sprintf(symlink_fname, "%s.symlink", fname);

        /* ensure symlink_fname does not exist */
        unlink(symlink_fname);

        fp = fopen(fname, "w");
        if (fp == NULL) {
            err = -1;
            ERR_HANDLE("fopen", err)
        }
        else if (fclose(fp) != 0) {
            err = -1;
            ERR_HANDLE("fclose", err)
        }
        else if (symlink(fname, symlink_fname) != 0) {
            /* create a symbolic link to the regular file */
            err = -1;
            ERR_HANDLE("symlink", err)
        }
        sync();
    }
    MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (err != 0) {
        nerrs++;
        goto fn_exit;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /* symlink_fname may have file system type prefix */
    sprintf(symlink_fname, "%s.symlink", filename);

    /* create a file in NC_CLOBBER mode */
    err = ncmpi_create(MPI_COMM_WORLD, symlink_fname, NC_CLOBBER, MPI_INFO_NULL, &ncid);
    CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        /* when calling POSIX I/O, remove file type prefix from file name */
        sprintf(symlink_fname, "%s.symlink", fname);

        /* check if symlink_fname is still a symbolic link */
        err = lstat(symlink_fname, &statbuf);
        ERR_HANDLE("lstat", err)

        if (S_ISLNK(statbuf.st_mode)) {
            if (verbose)
                printf("file \"%s\" is a symbolic link\n", symlink_fname);
        }
        else {
            nerrs++;
            fprintf(stderr,"file \"%s\" is NOT a symbolic link\n", symlink_fname);
        }

        /* delete test files */
        err = unlink(symlink_fname);
        ERR_HANDLE("unlink", err)
        sync();
    }

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
    free(symlink_fname);

fn_exit:
    MPI_Finalize();
    return (nerrs > 0);
}

