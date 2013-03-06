/*********************************************************************
 *
 *  Copyright (C) 2013, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <pnetcdf.h>

#define NVARS 8
#define NX 5

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests if the file header size and variable offsets are properly
 * set when using a different set of alignment hints to open an existing file
 * and entering the redef mode to add more dimensions, attributes, and
 * variables.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o alignment_test alignment_test.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 alignment_test testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define ERR {if(err!=NC_NOERR)printf("Error at line=%d: %s\n", __LINE__, ncmpi_strerror(err));}

int main(int argc, char** argv) {
    int i, j, rank, nprocs, err, verbose=0, nfailed=0, nfailed_all;
    int ncid, cmode, varid[NVARS], dimid[2], *buf;
    char str[32];
    MPI_Offset start[2], count[2];
    MPI_Offset new_var_off[NVARS*2], old_var_off[NVARS*2];
    MPI_Offset header_size[2], header_extent[2];
    MPI_Info info=MPI_INFO_NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc != 2) {
        if (!rank) printf("Usage: %s filename\n",argv[0]);
        MPI_Finalize();
        return 0;
    }

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, argv[1], cmode, info, &ncid); ERR

    /* define dimension */
    err = ncmpi_def_dim(ncid, "Y", NC_UNLIMITED, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "X", NX*nprocs, &dimid[1]); ERR

    /* Odd numbers are fixed variables, even numbers are record variables */
    for (i=0; i<NVARS; i++) {
        if (i%2) {
            sprintf(str,"fixed_var_%d",i);
            err = ncmpi_def_var(ncid, str, NC_INT, 1, dimid+1, &varid[i]); ERR
        }
        else {
            sprintf(str,"record_var_%d",i);
            err = ncmpi_def_var(ncid, str, NC_INT, 2, dimid, &varid[i]); ERR
        }
    }
    err = ncmpi_enddef(ncid); ERR

    /* write all variables */
    buf = (int*) malloc(NX * sizeof(int));
    for (i=0; i<NVARS; i++) {
        for (j=0; j<NX; j++) buf[j] = rank*1000 + i*10 + j;
        if (i%2) {
            start[0] = NX*rank;
            count[0] = NX;
            err = ncmpi_put_vara_int_all(ncid, varid[i], start, count, buf); ERR
        }
        else {
            start[0] = 0; start[1] = NX*rank;
            count[0] = 1; count[1] = NX;
            err = ncmpi_put_vara_int_all(ncid, varid[i], start, count, buf); ERR
            for (j=0; j<NX; j++) buf[j] = rank*1000 + 100 + i*10 + j;
            start[0] = 1; /* write 2nd record */
            err = ncmpi_put_vara_int_all(ncid, varid[i], start, count, buf); ERR
        }
    }
    err = ncmpi_close(ncid); ERR

    /* Now, reopen the file and grow the header and read data back */

    /* mimic netCDF that does not do alignments */
    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_header_align_size", "1"); /* size in bytes */
    MPI_Info_set(info, "nc_var_align_size",    "1"); /* size in bytes */

    /* open the file for adding more metadata */
    err = ncmpi_open(MPI_COMM_WORLD, argv[1], NC_WRITE, info, &ncid); ERR

    /* get header size and extent, and offsets of all variables */
    err = ncmpi_inq_header_size(ncid, &header_size[0]); ERR
    err = ncmpi_inq_header_extent(ncid, &header_extent[0]); ERR
    for (i=0; i<NVARS; i++) {
        err = ncmpi_inq_varoffset(ncid, varid[i], &old_var_off[i]); ERR
    }

    /* enter redef mode */
    err = ncmpi_redef(ncid); ERR

    /* add attributes to make header grow */
    for (i=0; i<NVARS; i++) {
        sprintf(str, "annotation_for_var_%d",i);
        err = ncmpi_put_att_text(ncid, varid[i], "text_attr", strlen(str), str); ERR
    }

    /* add new dimensions */
    int new_dimid[3];
    err = ncmpi_def_dim(ncid, "new_dim_a", 5,         &new_dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "new_dim_b", 4,         &new_dimid[1]); ERR
    err = ncmpi_def_dim(ncid, "new_dim_c", NX*nprocs, &new_dimid[2]); ERR

    /* add new variables */
    int new_varid[NVARS];
    for (i=0; i<NVARS; i++) {
        if (i%2) {
            sprintf(str,"record_var_%d",i+NVARS);
            err = ncmpi_def_var(ncid, str, NC_INT, 2, dimid, &new_varid[i]); ERR
        }
        else {
            sprintf(str,"fixed_var_%d",i+NVARS);
            err = ncmpi_def_var(ncid, str, NC_INT, 1, new_dimid+2, &new_varid[i]); ERR
        }
    }
    err = ncmpi_enddef(ncid); ERR

    /* get the new header size and extent, also all variables' starting
       file offsets */
    err = ncmpi_inq_header_size(ncid, &header_size[1]); ERR
    err = ncmpi_inq_header_extent(ncid, &header_extent[1]); ERR
    if (rank == 0 && verbose) {
        printf("NX = %d (integer type)\n",NX);
        printf("old header_size  =%lld new header_size  =%lld\n",header_size[0],header_size[1]);
        printf("old header_extent=%lld new header_extent=%lld\n",header_extent[0],header_extent[1]);
        for (i=0; i<2*NVARS; i++) {
            int vid = (i < NVARS) ? varid[i] : new_varid[i-NVARS];
            err = ncmpi_inq_varoffset(ncid, vid, &new_var_off[i]); ERR
            if (i < NVARS)
                printf("var[%2d] old offset=%4lld new offset=%4lld\n",i,old_var_off[i],new_var_off[i]);
            else
                printf("var[%2d] additional new offset=%4lld\n",i,new_var_off[i]);
        }
    }

    /* write to the new variables */
    for (i=0; i<NVARS; i++) {
        for (j=0; j<NX; j++) buf[j] = -1 * (i*10 + j);
        if (i%2) {
            start[0] = 0; start[1] = NX*rank;
            count[0] = 1; count[1] = NX;
            err = ncmpi_put_vara_int_all(ncid, new_varid[i], start, count, buf); ERR
            for (j=0; j<NX; j++) buf[j] = -1 * (100 + i*10 + j);
            start[0] = 1; /* write 2nd record */
            err = ncmpi_put_vara_int_all(ncid, new_varid[i], start, count, buf); ERR
        }
        else {
            start[0] = NX*rank;
            count[0] = NX;
            err = ncmpi_put_vara_int_all(ncid, new_varid[i], start, count, buf); ERR
        }
    }

    /* read old variables and check their contents */
    for (i=0; i<NVARS; i++) {
        if (i%2) {
            start[0] = NX*rank;
            count[0] = NX;
            err = ncmpi_get_vara_int_all(ncid, varid[i], start, count, buf); ERR
            for (j=0; j<NX; j++)
                if (buf[j] != rank*1000 + i*10 + j) {
                    printf("read error i=%d buf[j=%d]=%d != %d\n",i,j,buf[j],rank*1000+i*10+j);
                    nfailed++;
                }
        }
        else {
            start[0] = 0; start[1] = NX*rank;
            count[0] = 1; count[1] = NX;
            err = ncmpi_get_vara_int_all(ncid, varid[i], start, count, buf); ERR
            for (j=0; j<NX; j++)
                if (buf[j] != rank*1000+i*10+j) {
                    printf("read error i=%d buf[j=%d]=%d != %d\n",i,j,buf[j],rank*1000+i*10+j);
                    nfailed++;
                }
            start[0] = 1;
            err = ncmpi_get_vara_int_all(ncid, varid[i], start, count, buf); ERR
            for (j=0; j<NX; j++)
                if (buf[j] != rank*1000 + 100 + i*10 + j) {
                    printf("read error i=%d buf[j=%d]=%d != %d\n",i,j,buf[j],rank*1000+100+i*10+j);
                    nfailed++;
                }
        }
    }
    err = ncmpi_close(ncid); ERR
    MPI_Info_free(&info);

    free(buf);
    MPI_Reduce(&nfailed, &nfailed_all, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0 && verbose) {
        if (nfailed_all == 0) printf("alignment test (%s) --- passed\n",argv[0]);
        else printf("alignment test (%s) failed with %d mismatches\n",argv[0],nfailed_all);
    }
    MPI_Finalize();
    return 0;
}

