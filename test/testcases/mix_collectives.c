/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    char filename[256];
    int i, j, err, nerrs=0, rank, nprocs;
    int ncid, dimid[2], varid, varids[4];
    MPI_Offset start[2], count[2], stride[2], imap[2];
    int   *check_buf, buf[6][4];
    int   g_buf[96] = {
    NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, 100, 101, 102, 103,
    NC_FILL_INT, 0,           NC_FILL_INT, NC_FILL_INT, 104, 105, 106, 107,
    NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, 108, 109, 110, 111,
    NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, 112, 113, 114, 115,
    NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, 116, 117, 118, 119,
    NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, 120, 121, 122, 123,
    200,         NC_FILL_INT,         201, NC_FILL_INT, 300, 306, 312, 318,
    NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, 301, 307, 313, 319,
    202,         NC_FILL_INT,         203, NC_FILL_INT, 302, 308, 314, 320,
    NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, 303, 309, 315, 321,
    204,         NC_FILL_INT,         205, NC_FILL_INT, 304, 310, 316, 322,
    NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, NC_FILL_INT, 305, 311, 317, 323};

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
        sprintf(cmd_str, "*** TESTING C   %s for get/put varm ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

#ifdef DEBUG
    if (nprocs > 1 && rank == 0)
        printf("Warning: %s is designed to run on 1 process\n", argv[0]);
#endif

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER | NC_64BIT_DATA,
                       MPI_INFO_NULL, &ncid); CHECK_ERR

    /* define a variable of a 6 x 4 integer array in the nc file */
    err = ncmpi_def_dim(ncid, "Y", 12, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 8, &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var", NC_INT, 2, dimid, &varid); CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    for (j=0; j<6; j++) for (i=0; i<4; i++) buf[j][i] = j*4+i + rank*100;

    /* now, each of 4 processes makes a call to different kinds of put APIs */
    if (rank == 0) {
        /* write a 1 x 1 subarray with contents:
               _, _, _, _, _, _, _, _,
               _, 0, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _
         */
        start[0] = 1; start[1] = 1;
        err = ncmpi_put_var1_int_all(ncid, varid, start, &buf[0][0]); CHECK_ERR
    }
    else if (rank == 1) {
        /* write a 6 x 4 subarray with contents:
               _, _, _, _, 100, 101, 102, 103,
               _, _, _, _, 104, 105, 106, 107,
               _, _, _, _, 108, 109, 110, 111,
               _, _, _, _, 112, 113, 114, 115,
               _, _, _, _, 116, 117, 118, 119,
               _, _, _, _, 120, 121, 122, 123,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _
         */
        start[0] = 0; start[1] = 4;
        count[0] = 6; count[1] = 4;
        err = ncmpi_put_vara_int_all(ncid, varid, start, count, &buf[0][0]); CHECK_ERR
    }
    else if (rank == 2) {
        /* write a strided 6 x 4 subarray with contents:
                 _, _,   _, _, _, _, _, _,
                 _, _,   _, _, _, _, _, _,
                 _, _,   _, _, _, _, _, _,
                 _, _,   _, _, _, _, _, _,
                 _, _,   _, _, _, _, _, _,
                 _, _,   _, _, _, _, _, _,
               200, _, 201, _, _, _, _, _,
                 _, _,   _, _, _, _, _, _,
               202, _, 203, _, _, _, _, _,
                 _, _,   _, _, _, _, _, _,
               204, _, 205, _, _, _, _, _,
                 _, _,   _, _, _, _, _, _
         */
        start[0]  = 6; start[1]  = 0;
        count[0]  = 3; count[1]  = 2;
        stride[0] = 2; stride[1] = 2;
        err = ncmpi_put_vars_int_all(ncid, varid, start, count, stride, &buf[0][0]); CHECK_ERR
    }
    else if (rank == 3) {
        /* write a 6 x 4 transported subarray with contents:
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, _, _, _, _,
               _, _, _, _, 300, 306, 312, 318,
               _, _, _, _, 301, 307, 313, 319,
               _, _, _, _, 302, 308, 314, 320,
               _, _, _, _, 303, 309, 315, 321,
               _, _, _, _, 304, 310, 316, 322,
               _, _, _, _, 305, 311, 317, 323

         */
        start[0]  = 6; start[1]  = 4;
        count[0]  = 6; count[1]  = 4;
        stride[0] = 1; stride[1] = 1;
        imap[0]   = 1; imap[1]   = 6;   /* would be {4, 1} if not transposing */
        err = ncmpi_put_varm_int_all(ncid, varid, start, count, stride, imap, &buf[0][0]); CHECK_ERR
    }
    else {
        start[0] = 6; start[1] = 4;
        count[0] = 0; count[1] = 0;
        err = ncmpi_put_vara_all(ncid, varid, start, count, &buf[0][0], 0, MPI_DATATYPE_NULL); CHECK_ERR
    }

    err = ncmpi_close(ncid); CHECK_ERR

    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_WRITE, MPI_INFO_NULL, &ncid); CHECK_ERR

    err = ncmpi_inq_varid(ncid, "var", &varid); CHECK_ERR

    check_buf = (int*) calloc(96, sizeof(int));

    err = ncmpi_get_var_int_all(ncid, varid, check_buf); CHECK_ERR

    /* read the whole variable and check contents */
    if (nprocs == 4) {
        for (i=0; i<96; i++) {
            if (check_buf[i] != g_buf[i]) {
#ifdef PRINT_ERR_ON_SCREEN
                printf("Error at line %d in %s: expecting var[%d]=%d but got %d\n",
                __LINE__,__FILE__,i,g_buf[i],check_buf[i]);
#endif
                nerrs++;
                break;
            }
        }
    }
    free(check_buf);

    /* now, each of 4 processes makes a call to different kinds of get APIs */
    for (j=0; j<6; j++) for (i=0; i<4; i++) buf[j][i] = -1;
    if (rank == 0) {
        start[0] = 1; start[1] = 1;
        err = ncmpi_get_var1_int_all(ncid, varid, start, &buf[0][0]); CHECK_ERR
        if (buf[0][0] != rank) {
#ifdef PRINT_ERR_ON_SCREEN
            printf("Error at line %d in %s: expecting buf[0][0]=%d but got %d\n",
            __LINE__,__FILE__,rank,buf[0][0]);
#endif
            nerrs++;
        }
    }
    else if (rank == 1) {
        start[0] = 0; start[1] = 4;
        count[0] = 6; count[1] = 4;
        err = ncmpi_get_vara_int_all(ncid, varid, start, count, &buf[0][0]); CHECK_ERR
        for (j=0; j<6; j++) for (i=0; i<4; i++) {
            if (buf[j][i] != j*4+i + rank*100) {
#ifdef PRINT_ERR_ON_SCREEN
                printf("Error at line %d in %s: expecting var[%d]=%d but got %d\n",
                __LINE__,__FILE__,i, j*4+i + rank*100, buf[j][i]);
#endif
                nerrs++;
            }
        }
    }
    else if (rank == 2) {
        start[0]  = 6; start[1]  = 0;
        count[0]  = 3; count[1]  = 2;
        stride[0] = 2; stride[1] = 2;
        err = ncmpi_get_vars_int_all(ncid, varid, start, count, stride, &buf[0][0]); CHECK_ERR
        int *val = &buf[0][0];
        for (j=0; j<count[0]*count[1]; j++) {
            if (*val != j + rank*100) {
#ifdef PRINT_ERR_ON_SCREEN
                printf("Error at line %d in %s: expecting var[%d]=%d but got %d\n",
                __LINE__,__FILE__,j, j+rank*100, *val);
#endif
                nerrs++;
            }
            val++;
        }
    }
    else if (rank == 3) {
        start[0]  = 6; start[1]  = 4;
        count[0]  = 6; count[1]  = 4;
        stride[0] = 1; stride[1] = 1;
        imap[0]   = 1; imap[1]   = 6;   /* would be {4, 1} if not transposing */
        err = ncmpi_get_varm_int_all(ncid, varid, start, count, stride, imap, &buf[0][0]); CHECK_ERR
        for (j=0; j<6; j++) for (i=0; i<4; i++) {
            if (buf[j][i] != j*4+i + rank*100) {
#ifdef PRINT_ERR_ON_SCREEN
                printf("Error at line %d in %s: expecting var[%d][%d]=%d but got %d\n",
                __LINE__,__FILE__,j,i, j*4+i + rank*100, buf[j][i]);
#endif
                nerrs++;
            }
        }
    }
    else {
        start[0] = 6; start[1] = 4;
        count[0] = 0; count[1] = 0;
        err = ncmpi_get_vara_all(ncid, varid, start, count, &buf[0][0], 0, MPI_DATATYPE_NULL); CHECK_ERR
        for (j=0; j<6; j++) for (i=0; i<4; i++) {
            if (buf[j][i] != -1) {
#ifdef PRINT_ERR_ON_SCREEN
                printf("Error at line %d in %s: expecting var[%d][%d]=%d but got %d\n",
                __LINE__,__FILE__,j,i, -1, buf[j][i]);
#endif
                nerrs++;
            }
        }
    }

    /* test when different processes call put APIs with different varid */
    err = ncmpi_redef(ncid); CHECK_ERR
    err = ncmpi_def_var(ncid, "scalar0", NC_INT, 0, NULL, &varids[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "scalar1", NC_INT, 0, NULL, &varids[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "scalar2", NC_INT, 0, NULL, &varids[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "scalar3", NC_INT, 0, NULL, &varids[3]); CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (rank < 4) {
        err = ncmpi_put_var_int_all(ncid, varids[rank], &rank); CHECK_ERR
    }
    else { /* make zero-length request */
        start[0] = 0; start[1] = 0;
        count[0] = 0; count[1] = 0;
        err = ncmpi_put_vara_int_all(ncid, varid, start, count, &rank); CHECK_ERR
    }

    err = ncmpi_close(ncid); CHECK_ERR

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

