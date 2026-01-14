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

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    char *symlink_fname, *fname;
    int  err, nerrs=0, rank, ncid, verbose=0;
    struct stat statbuf;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* remove file system type prefix substring */
    fname = remove_file_system_type_prefix(out_path);

    symlink_fname = (char*) malloc(strlen(out_path) + 10);

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
    if (err != 0) return 1;
    MPI_Barrier(MPI_COMM_WORLD);

    /* symlink_fname may have file system type prefix */
    sprintf(symlink_fname, "%s.symlink", out_path);

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a file in NC_CLOBBER mode */
    err = ncmpi_create(MPI_COMM_WORLD, symlink_fname, NC_CLOBBER, info, &ncid);
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

    free(symlink_fname);

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_64BIT_DATA};
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 0; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 0; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "NC_CLOBBER on symlink file", opt, test_io);

    MPI_Finalize();

    return err;
}
