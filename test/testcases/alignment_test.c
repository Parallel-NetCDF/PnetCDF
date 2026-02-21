/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests if the file header size and variable offsets are properly
 * set when using a different set of alignment hints to open an existing file
 * and entering the redef mode to add more dimensions, attributes, and
 * variables, causing the expansion of the header.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o alignment_test alignment_test.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 alignment_test testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <assert.h>
#include <mpi.h>
#include <pnetcdf.h>

#include <testutils.h>

#define NVARS 8
#define NX 70

#define TEST_FIXED_VAR
#define TEST_RECORD_VAR

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    global_info)
{
    int i, j, rank, nprocs, err, verbose=0, nerrs=0;
    int ncid, varid[NVARS], dimid[2], *buf;
    char str[32];
    MPI_Offset start[2], count[2];
    MPI_Offset new_var_off[NVARS*2], old_var_off[NVARS*2];
    MPI_Offset header_size[2], header_extent[2];
    MPI_Info info=MPI_INFO_NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    MPI_Info_dup(global_info, &info);

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file for writing ----------------------------------------*/
    err = ncmpi_create(MPI_COMM_WORLD, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERROUT

    /* define dimension */
    err = ncmpi_def_dim(ncid, "Y", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", NX*nprocs, &dimid[1]); CHECK_ERR

    /* Odd numbers are fixed variables, even numbers are record variables */
    for (i=0; i<NVARS; i++) {
#ifdef TEST_FIXED_VAR
        if (i%2) {
            sprintf(str,"fixed_var_%d",i);
            err = ncmpi_def_var(ncid, str, NC_INT, 1, dimid+1, &varid[i]); CHECK_ERR
        }
#endif
#ifdef TEST_RECORD_VAR
        if (i%2 == 0) {
            sprintf(str,"record_var_%d",i);
            err = ncmpi_def_var(ncid, str, NC_INT, 2, dimid, &varid[i]); CHECK_ERR
        }
#endif
    }
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* write all variables */
    buf = (int*) malloc(sizeof(int) * NX);
    for (i=0; i<NVARS; i++) {
        for (j=0; j<NX; j++) buf[j] = rank*1000 + i*10 + j;
#ifdef TEST_FIXED_VAR
        if (i%2) {
            start[0] = NX*rank;
            count[0] = NX;
            if (coll_io)
                err = ncmpi_put_vara_int_all(ncid, varid[i], start, count, buf);
            else
                err = ncmpi_put_vara_int(ncid, varid[i], start, count, buf);
            CHECK_ERR

            /* check if user put buffer contents altered */
            for (j=0; j<NX; j++) {
                if (buf[j] != rank*1000 + i*10 + j) {
                    printf("Error at line %d in %s: user put buffer[%d] altered from %d to %d\n",
                           __LINE__,__FILE__,j, rank*1000 + i*10 + j, buf[j]);
                    nerrs++;
                }
            }
        }
#endif
#ifdef TEST_RECORD_VAR
        if (i%2 == 0) {
            start[0] = 0; start[1] = NX*rank;
            count[0] = 1; count[1] = NX;
            if (coll_io)
                err = ncmpi_put_vara_int_all(ncid, varid[i], start, count, buf);
            else
                err = ncmpi_put_vara_int(ncid, varid[i], start, count, buf);
            CHECK_ERR
            for (j=0; j<NX; j++) buf[j] = rank*1000 + 100 + i*10 + j;
            start[0] = 1; /* write 2nd record */
            if (coll_io)
                err = ncmpi_put_vara_int_all(ncid, varid[i], start, count, buf);
            else
                err = ncmpi_put_vara_int(ncid, varid[i], start, count, buf);
            CHECK_ERR
            /* check if user put buffer contents altered */
            for (j=0; j<NX; j++) {
                if (buf[j] != rank*1000 + 100 + i*10 + j) {
                    printf("Error at line %d in %s: user put buffer[%d] altered from %d to %d\n",
                           __LINE__,__FILE__,j, rank*1000 + 100 + i*10 + j, buf[j]);
                    nerrs++;
                }
            }
        }
#endif
    }

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    err = ncmpi_close(ncid); CHECK_ERR

    /* Now, reopen the file and grow the header and read data back */

    /* mimic netCDF that does not do alignments */
    MPI_Info_set(info, "nc_var_align_size", "197"); /* size in bytes */

    /* open the file for adding more metadata */
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_WRITE, info, &ncid);
    CHECK_ERROUT

    /* get header size and extent */
    err = ncmpi_inq_header_size(ncid, &header_size[0]); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &header_extent[0]); CHECK_ERR
    /* get offsets of all variables */
    for (i=0; i<NVARS; i++) {
#ifdef TEST_FIXED_VAR
        if (i%2)
            err = ncmpi_inq_varoffset(ncid, varid[i], &old_var_off[i]);
#endif
#ifdef TEST_RECORD_VAR
        if (i%2==0)
            err = ncmpi_inq_varoffset(ncid, varid[i], &old_var_off[i]);
#endif
        CHECK_ERR
    }

    /* re-enter define mode */
    err = ncmpi_redef(ncid); CHECK_ERR

    /* add attributes to make header grow */
    for (i=0; i<NVARS; i++) {
        sprintf(str, "annotation_for_var_%d",i);
#ifdef TEST_FIXED_VAR
        if (i%2)
            err = ncmpi_put_att_text(ncid, varid[i], "text_attr", strlen(str), str);
#endif
#ifdef TEST_RECORD_VAR
        if (i%2==0)
            err = ncmpi_put_att_text(ncid, varid[i], "text_attr", strlen(str), str);
#endif
        CHECK_ERR
    }

    /* add new dimensions */
    int new_dimid[3];
    err = ncmpi_def_dim(ncid, "new_dim_a", 5,         &new_dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "new_dim_b", 4,         &new_dimid[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "new_dim_c", NX*nprocs, &new_dimid[2]); CHECK_ERR

    /* add new variables */
    int new_varid[NVARS];
    for (i=0; i<NVARS; i++) {
#ifdef TEST_FIXED_VAR
        if (i%2 == 0) {
            sprintf(str,"fixed_var_%d",i+NVARS);
            err = ncmpi_def_var(ncid, str, NC_INT, 1, new_dimid+2, &new_varid[i]); CHECK_ERR
        }
#endif
#ifdef TEST_RECORD_VAR
        if (i%2 == 1) {
            sprintf(str,"record_var_%d",i+NVARS);
            err = ncmpi_def_var(ncid, str, NC_INT, 2, dimid, &new_varid[i]); CHECK_ERR
        }
#endif
    }
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* get the new header size and extent, also all variables' starting
       file offsets */
    err = ncmpi_inq_header_size(ncid, &header_size[1]); CHECK_ERR
    err = ncmpi_inq_header_extent(ncid, &header_extent[1]); CHECK_ERR
    if (rank == 0 && verbose) {
        printf("NX = %d (integer type)\n",NX);
        printf("old header_size  ="OFFFMT" new header_size  ="OFFFMT"\n",header_size[0],header_size[1]);
        printf("old header_extent="OFFFMT" new header_extent="OFFFMT"\n",header_extent[0],header_extent[1]);

#ifdef TEST_FIXED_VAR
        for (i=1; i<NVARS; i+=2) {
            err = ncmpi_inq_varoffset(ncid, varid[i], &new_var_off[i]); CHECK_ERR
            printf("old fixed  var[%2d] old offset="OFFFMT" new offset="OFFFMT"\n",i,old_var_off[i],new_var_off[i]);
        }
        for (i=NVARS; i<2*NVARS; i++) {
            if (i%2 == 0) {
                err = ncmpi_inq_varoffset(ncid, new_varid[i-NVARS], &new_var_off[i]); CHECK_ERR
                printf("new fixed  var[%2d]                 new offset="OFFFMT"\n",i,new_var_off[i]);
            }
        }
#endif
#ifdef TEST_RECORD_VAR
        for (i=0; i<NVARS; i+=2) {
            err = ncmpi_inq_varoffset(ncid, varid[i], &new_var_off[i]); CHECK_ERR
            printf("old record var[%2d] old offset="OFFFMT" new offset="OFFFMT"\n",i,old_var_off[i],new_var_off[i]);
        }
        for (i=NVARS; i<2*NVARS; i++) {
            if (i%2) {
                err = ncmpi_inq_varoffset(ncid, new_varid[i-NVARS], &new_var_off[i]); CHECK_ERR
                printf("new record var[%2d]                 new offset="OFFFMT"\n",i,new_var_off[i]);
            }
        }
#endif
    }

    /* write to the new variables */
    for (i=0; i<NVARS; i++) {
        for (j=0; j<NX; j++) buf[j] = -1 * (i*10 + j);
#ifdef TEST_FIXED_VAR
        if (i%2 == 0) {
            start[0] = NX*rank;
            count[0] = NX;
            if (coll_io)
                err = ncmpi_put_vara_int_all(ncid, new_varid[i], start, count, buf);
            else
                err = ncmpi_put_vara_int(ncid, new_varid[i], start, count, buf);
            CHECK_ERR
            sprintf(str,"fixed_var_%d",i);
            /* check if user put buffer contents altered */
            for (j=0; j<NX; j++) {
                if (buf[j] != -1 * (i*10 + j)) {
                    printf("Error at %d in %s: var %s put buffer[%d] altered from %d to %d\n",
                           __LINE__,__FILE__, str, j, -1 * (i*10 + j), buf[j]);
                    nerrs++;
                }
            }
        }
#endif
#ifdef TEST_RECORD_VAR
        if (i%2 == 1) {
            start[0] = 0; start[1] = NX*rank;
            count[0] = 1; count[1] = NX;
            if (coll_io)
                err = ncmpi_put_vara_int_all(ncid, new_varid[i], start, count, buf);
            else
                err = ncmpi_put_vara_int(ncid, new_varid[i], start, count, buf);
            CHECK_ERR
            for (j=0; j<NX; j++) buf[j] = -1 * (100 + i*10 + j);
            start[0] = 1; /* write 2nd record */
            if (coll_io)
                err = ncmpi_put_vara_int_all(ncid, new_varid[i], start, count, buf);
            else
                err = ncmpi_put_vara_int(ncid, new_varid[i], start, count, buf);
            CHECK_ERR
            sprintf(str,"record_var_%d",i);
            /* check if user put buffer contents altered */
            for (j=0; j<NX; j++) {
                if (buf[j] != -1 * (100 + i*10 + j)) {
                    printf("Error at %d in %s: var %s put buffer[%d] altered from %d to %d\n",
                           __LINE__,__FILE__, str, j, -1 * (100 + i*10 + j), buf[j]);
                    nerrs++;
                }
            }
        }
#endif
    }

    /* file sync before reading */
    err = ncmpi_sync(ncid);
    CHECK_ERR
    MPI_Barrier(MPI_COMM_WORLD);

    /* read old variables and check their contents */
    for (i=0; i<NVARS; i++) {
#ifdef TEST_FIXED_VAR
        if (i%2) {
            start[0] = NX*rank;
            count[0] = NX;
            for (j=0; j<NX; j++) buf[j] = -1;
            if (coll_io)
                err = ncmpi_get_vara_int_all(ncid, varid[i], start, count, buf);
            else
                err = ncmpi_get_vara_int(ncid, varid[i], start, count, buf);
            CHECK_ERR
            sprintf(str,"fixed_var_%d",i);
            for (j=0; j<NX; j++)
                if (buf[j] != rank*1000 + i*10 + j) {
                    printf("Error at %d: var %s i=%d buf[j=%d]=%d != %d\n",
                           __LINE__,str,i,j,buf[j],rank*1000+i*10+j);
                    nerrs++;
                    break;
                }
        }
#endif
#ifdef TEST_RECORD_VAR
        if (i%2 == 0) {
            start[0] = 0; start[1] = NX*rank;
            count[0] = 1; count[1] = NX;
            if (coll_io)
                err = ncmpi_get_vara_int_all(ncid, varid[i], start, count, buf);
            else
                err = ncmpi_get_vara_int(ncid, varid[i], start, count, buf);
            CHECK_ERR
            for (j=0; j<NX; j++)
                if (buf[j] != rank*1000+i*10+j) {
                    printf("read error at %d: i=%d buf[j=%d]=%d != %d\n",__LINE__,i,j,buf[j],rank*1000+i*10+j);
                    nerrs++;
                    break;
                }
            start[0] = 1;
            if (coll_io)
                err = ncmpi_get_vara_int_all(ncid, varid[i], start, count, buf);
            else
                err = ncmpi_get_vara_int(ncid, varid[i], start, count, buf);
            CHECK_ERR
            sprintf(str,"record_var_%d",i);
            for (j=0; j<NX; j++)
                if (buf[j] != rank*1000 + 100 + i*10 + j) {
                    printf("Error at %d: var %s i=%d buf[j=%d]=%d != %d\n",
                           __LINE__, str,i,j,buf[j],rank*1000+100+i*10+j);
                    nerrs++;
                    break;
                }
        }
#endif
    }
    err = ncmpi_close(ncid); CHECK_ERR
    free(buf);

    if (info != MPI_INFO_NULL) MPI_Info_free(&info);

err_out:
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
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 1; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "alignment hints", opt, test_io);

    MPI_Finalize();

    return err;
}
