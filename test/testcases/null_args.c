/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* This program tests whether the correct error codes can be returned when
 * using NULL arguments for start, count, stride, or imap
 */

#include <stdio.h>
#include <stdlib.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define EXP_ERR_MSG(exp,msg) { \
    if (err != exp) { \
        nerrs++; \
        fprintf(stderr, "Error at line %d in %s: (%s) expect %s but got %s\n", \
                __LINE__, __FILE__, msg, \
                ncmpi_strerrno(exp), ncmpi_strerrno(err)); \
    } \
}

int main(int argc, char **argv)
{
    char filename[256];
    int err, nerrs=0, ncid, dimid[2], varid;
    int nprocs, rank, buf[100];
    MPI_Offset start[2], count[2], stride[2], imap[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for NULL arguments ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    err = ncmpi_create(MPI_COMM_WORLD, filename, 0, MPI_INFO_NULL, &ncid);
    EXP_ERR_MSG(NC_NOERR, "create")
    err = ncmpi_def_dim(ncid, "Y", NC_UNLIMITED, &dimid[0]);
    EXP_ERR_MSG(NC_NOERR,"def_dim Y")
    err = ncmpi_def_dim(ncid, "X", 10, &dimid[1]);
    EXP_ERR_MSG(NC_NOERR,"def_dim X")
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimid, &varid);
    EXP_ERR_MSG(NC_NOERR,"def_var")
    err = ncmpi_enddef(ncid);
    EXP_ERR_MSG(NC_NOERR,"enddef")
    err = ncmpi_begin_indep_data(ncid);
    EXP_ERR_MSG(NC_NOERR,"begin_indep_data")

     start[0] =  start[1] = 0;
     count[0] =  count[1] = 1;
    stride[0] = stride[1] = 1;
      imap[0] =   imap[1] = 1;

    memset(buf, 0, 100*sizeof(int));

    /*---- test put_var1 ---- */
    err = ncmpi_put_var1_int(ncid, varid, start, buf);
    EXP_ERR_MSG(NC_NOERR, "put_var1")

    err = ncmpi_put_var1_int(ncid, varid, NULL, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "put_var1 start=NULL")

    /*---- test put_vara ---- */
    err = ncmpi_put_vara_int(ncid, varid, start, count, buf);
    EXP_ERR_MSG(NC_NOERR, "put_vara")

    err = ncmpi_put_vara_int(ncid, varid, NULL, count, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "put_vara start=NULL")

    err = ncmpi_put_vara_int(ncid, varid, start, NULL, buf);
    EXP_ERR_MSG(NC_EEDGE, "put_vara count=NULL")

    err = ncmpi_put_vara_int(ncid, varid, NULL, NULL, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "put_vara start=count=NULL")

    /*---- test put_vars ---- */
    err = ncmpi_put_vars_int(ncid, varid, start, count, stride, buf);
    EXP_ERR_MSG(NC_NOERR, "put_vars")

    err = ncmpi_put_vars_int(ncid, varid, NULL, count, stride, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "put_vars start=NULL")

    err = ncmpi_put_vars_int(ncid, varid, start, NULL, stride, buf);
    EXP_ERR_MSG(NC_EEDGE, "put_vars count=NULL")

    err = ncmpi_put_vars_int(ncid, varid, start, count, NULL, buf);
    EXP_ERR_MSG(NC_NOERR, "put_vars stride=NULL")

    err = ncmpi_put_vars_int(ncid, varid, NULL, NULL, stride, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "put_vars start=count=NULL")

    err = ncmpi_put_vars_int(ncid, varid, NULL, count, NULL, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "put_vars start=stride=NULL")

    err = ncmpi_put_vars_int(ncid, varid, start, NULL, NULL, buf);
    EXP_ERR_MSG(NC_EEDGE, "put_vars count=stride=NULL")

    err = ncmpi_put_vars_int(ncid, varid, NULL, NULL, NULL, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "put_vars start=count=stride=NULL")

    /*---- test put_varm ---- */
    err = ncmpi_put_varm_int(ncid, varid, start, count, stride, imap, buf);
    EXP_ERR_MSG(NC_NOERR, "put_varm")

    err = ncmpi_put_varm_int(ncid, varid, NULL, count, stride, imap, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "put_varm start=NULL")

    err = ncmpi_put_varm_int(ncid, varid, start, NULL, stride, imap, buf);
    EXP_ERR_MSG(NC_EEDGE, "put_varm count=NULL")

    err = ncmpi_put_varm_int(ncid, varid, start, count, NULL, imap, buf);
    EXP_ERR_MSG(NC_NOERR, "put_varm stride=NULL")

    err = ncmpi_put_varm_int(ncid, varid, start, count, stride, NULL, buf);
    EXP_ERR_MSG(NC_NOERR, "put_varm imap=NULL")

    err = ncmpi_put_varm_int(ncid, varid, NULL, NULL, stride, imap, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "put_varm start=count=NULL")

    err = ncmpi_put_varm_int(ncid, varid, NULL, count, NULL, imap, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "put_varm start=stride=NULL")

    err = ncmpi_put_varm_int(ncid, varid, NULL, count, stride, NULL, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "put_varm start=imap=NULL")

    err = ncmpi_put_varm_int(ncid, varid, start, NULL, NULL, imap, buf);
    EXP_ERR_MSG(NC_EEDGE, "put_varm count=stride=NULL")

    err = ncmpi_put_varm_int(ncid, varid, start, NULL, stride, NULL, buf);
    EXP_ERR_MSG(NC_EEDGE, "put_varm count=imap=NULL")

    err = ncmpi_put_varm_int(ncid, varid, start, count, NULL, NULL, buf);
    EXP_ERR_MSG(NC_NOERR, "put_varm stride=imap=NULL")

    err = ncmpi_put_varm_int(ncid, varid, NULL, NULL, NULL, imap, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "put_varm start=count=stride=NULL")

    err = ncmpi_put_varm_int(ncid, varid, NULL, NULL, stride, NULL, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "put_varm start=count=imap=NULL")

    err = ncmpi_put_varm_int(ncid, varid, start, NULL, NULL, NULL, buf);
    EXP_ERR_MSG(NC_EEDGE, "put_varm count=stride=imap=NULL")

    err = ncmpi_put_varm_int(ncid, varid, NULL, NULL, NULL, NULL, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "put_varm start=count=stride=imap=NULL")

    /*---- test get_var1 ---- */
    err = ncmpi_get_var1_int(ncid, varid, start, buf);
    EXP_ERR_MSG(NC_NOERR, "get_var1")

    err = ncmpi_get_var1_int(ncid, varid, NULL, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "get_var1 start=NULL")

    /*---- test get_vara ---- */
    err = ncmpi_get_vara_int(ncid, varid, start, count, buf);
    EXP_ERR_MSG(NC_NOERR, "get_vara")

    err = ncmpi_get_vara_int(ncid, varid, NULL, count, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "get_vara start=NULL")

    err = ncmpi_get_vara_int(ncid, varid, start, NULL, buf);
    EXP_ERR_MSG(NC_EEDGE, "get_vara count=NULL")

    err = ncmpi_get_vara_int(ncid, varid, NULL, NULL, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "get_vara start=count=NULL")

    /*---- test get_vars ---- */
    err = ncmpi_get_vars_int(ncid, varid, start, count, stride, buf);
    EXP_ERR_MSG(NC_NOERR, "get_vars")

    err = ncmpi_get_vars_int(ncid, varid, NULL, count, stride, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "get_vars start=NULL")

    err = ncmpi_get_vars_int(ncid, varid, start, NULL, stride, buf);
    EXP_ERR_MSG(NC_EEDGE, "get_vars count=NULL")

    err = ncmpi_get_vars_int(ncid, varid, start, count, NULL, buf);
    EXP_ERR_MSG(NC_NOERR, "get_vars stride=NULL")

    err = ncmpi_get_vars_int(ncid, varid, NULL, NULL, stride, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "get_vars start=count=NULL")

    err = ncmpi_get_vars_int(ncid, varid, NULL, count, NULL, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "get_vars start=stride=NULL")

    err = ncmpi_get_vars_int(ncid, varid, start, NULL, NULL, buf);
    EXP_ERR_MSG(NC_EEDGE, "get_vars count=stride=NULL")

    err = ncmpi_get_vars_int(ncid, varid, NULL, NULL, NULL, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "get_vars start=count=stride=NULL")

    /*---- test get_varm ---- */
    err = ncmpi_get_varm_int(ncid, varid, start, count, stride, imap, buf);
    EXP_ERR_MSG(NC_NOERR, "get_varm")

    err = ncmpi_get_varm_int(ncid, varid, NULL, count, stride, imap, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "get_varm start=NULL")

    err = ncmpi_get_varm_int(ncid, varid, start, NULL, stride, imap, buf);
    EXP_ERR_MSG(NC_EEDGE, "get_varm count=NULL")

    err = ncmpi_get_varm_int(ncid, varid, start, count, NULL, imap, buf);
    EXP_ERR_MSG(NC_NOERR, "get_varm stride=NULL")

    err = ncmpi_get_varm_int(ncid, varid, start, count, stride, NULL, buf);
    EXP_ERR_MSG(NC_NOERR, "get_varm imap=NULL")

    err = ncmpi_get_varm_int(ncid, varid, NULL, NULL, stride, imap, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "get_varm start=count=NULL")

    err = ncmpi_get_varm_int(ncid, varid, NULL, count, NULL, imap, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "get_varm start=stride=NULL")

    err = ncmpi_get_varm_int(ncid, varid, NULL, count, stride, NULL, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "get_varm start=imap=NULL")

    err = ncmpi_get_varm_int(ncid, varid, start, NULL, NULL, imap, buf);
    EXP_ERR_MSG(NC_EEDGE, "get_varm count=stride=NULL")

    err = ncmpi_get_varm_int(ncid, varid, start, NULL, stride, NULL, buf);
    EXP_ERR_MSG(NC_EEDGE, "get_varm count=imap=NULL")

    err = ncmpi_get_varm_int(ncid, varid, start, count, NULL, NULL, buf);
    EXP_ERR_MSG(NC_NOERR, "get_varm stride=imap=NULL")

    err = ncmpi_get_varm_int(ncid, varid, NULL, NULL, NULL, imap, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "get_varm start=count=stride=NULL")

    err = ncmpi_get_varm_int(ncid, varid, NULL, NULL, stride, NULL, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "get_varm start=count=imap=NULL")

    err = ncmpi_get_varm_int(ncid, varid, start, NULL, NULL, NULL, buf);
    EXP_ERR_MSG(NC_EEDGE, "get_varm count=stride=imap=NULL")

    err = ncmpi_get_varm_int(ncid, varid, NULL, NULL, NULL, NULL, buf);
    EXP_ERR_MSG(NC_EINVALCOORDS, "get_varm start=count=stride=imap=NULL")

    err = ncmpi_close(ncid);
    EXP_ERR_MSG(NC_NOERR, "close")

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
