/*********************************************************************
 *
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program is to test
 *
 * a single large global attribute of size > NC_MAX_INT
 * a single large local  attribute of size > NC_MAX_INT
 *
 * Note the number of attributes is limited to NC_MAX_INT, because
 *    int ncmpi_inq_attid(int ncid, int varid, const char *name, int *idp);
 * argument idp is of type 4-byte int.
 * If this API is updated in the future to use large data type, then this test
 * can add more tests for number of attributes > NC_MAX_INT.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

int main(int argc, char** argv)
{
    char filename[256], *name, *buf;
    size_t i;
    int rank, nprocs, err, nerrs=0, verbose=0;
    int ncid, cmode, varid, dimid;
    MPI_Offset nelems, inq_nelems;
    MPI_Info info=MPI_INFO_NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* get command-line arguments */
    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) snprintf(filename, 256, "%s", argv[1]);
    else           strcpy(filename, "testfile.nc");
    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for one large ATTR", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    nelems = (MPI_Offset)NC_MAX_INT + 17;
    buf = (char*) malloc(nelems);

    /* create a new file and put a large global attribute -------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    /* put a large (> 2GiB) global attribute */
    for (i=0; i<nelems; i++) buf[i] = 'a' + i % 16;
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "large_attr", nelems, buf);
    if (!(cmode & NC_64BIT_DATA)) EXP_ERR(NC_EINVAL)
    else CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    if (nerrs > 0) goto err_out;

    /* open the file and read back the large global attribute ---------------*/
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, info, &ncid);
    CHECK_ERR
    if (err != NC_NOERR) goto err_out;

    err = ncmpi_inq_attlen(ncid, NC_GLOBAL, "large_attr", &inq_nelems);
    if (inq_nelems != nelems) {
        printf("Error at %s line %d: expecting attr nelems %lld but got %lld\n",
               __FILE__,__LINE__,nelems,inq_nelems);
        nerrs++;
    }

    for (i=0; i<nelems; i++) buf[i] = 0;
    err = ncmpi_get_att_text(ncid, NC_GLOBAL, "large_attr", buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) {
        char expect = 'a' + i % 16;
        if (buf[i] != 'a' + i % 16) {
            printf("Error at %s line %d: expecting attr[%zd] value %c but got %c\n",
                   __FILE__,__LINE__,i,expect,buf[i]);
            nerrs++;
            break;
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR
    if (nerrs > 0) goto err_out;

    /* create a new file and put a large local attribute -------------------*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid);
    CHECK_ERR

    err = ncmpi_def_var(ncid, "var", NC_INT, 1, &dimid, &varid);
    CHECK_ERR

    /* put a large (> 2GiB) global attribute */
    for (i=0; i<nelems; i++) buf[i] = 'a' + i % 16;
    err = ncmpi_put_att_text(ncid, varid, "large_attr", nelems, buf);
    if (!(cmode & NC_64BIT_DATA)) EXP_ERR(NC_EINVAL)
    else CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    if (nerrs > 0) goto err_out;

    /* open the file and read back the large global attribute ---------------*/
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, info, &ncid);
    CHECK_ERR
    if (err != NC_NOERR) goto err_out;

    err = ncmpi_inq_varid(ncid, "var", &varid);
    CHECK_ERR

    err = ncmpi_inq_attlen(ncid, varid, "large_attr", &inq_nelems);
    if (inq_nelems != nelems) {
        printf("Error at %s line %d: expecting attr len %lld but got %lld\n",
               __FILE__,__LINE__,nelems,inq_nelems);
        nerrs++;
    }

    for (i=0; i<nelems; i++) buf[i] = 0;
    err = ncmpi_get_att_text(ncid, varid, "large_attr", buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) {
        char expect = 'a' + i % 16;
        if (buf[i] != 'a' + i % 16) {
            printf("Error at %s line %d: expecting attr[%zd] value %c but got %c\n",
                   __FILE__,__LINE__,i,expect,buf[i]);
            nerrs++;
            break;
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR
    if (nerrs > 0) goto err_out;

    /* create a new file and put 2 global attributes, total size > 2 GiB ----*/
    nelems /= 2;

    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    /* put two global attributes (total size > 2GiB) */
    for (i=0; i<nelems; i++) buf[i] = 'a' + i % 16;
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "large_attr_0", nelems, buf);
    if (!(cmode & NC_64BIT_DATA)) EXP_ERR(NC_EINVAL)
    else CHECK_ERR

    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "large_attr_1", nelems, buf);
    if (!(cmode & NC_64BIT_DATA)) EXP_ERR(NC_EINVAL)
    else CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    if (nerrs > 0) goto err_out;

    /* open the file and read back the large global attributes --------------*/
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, info, &ncid);
    CHECK_ERR
    if (err != NC_NOERR) goto err_out;

    name = "large_attr_0";
    err = ncmpi_inq_attlen(ncid, NC_GLOBAL, name, &inq_nelems);
    if (inq_nelems != nelems) {
        printf("Error at %s line %d: expecting attr %s nelems %lld but got %lld\n",
               __FILE__,__LINE__,name, nelems,inq_nelems);
        nerrs++;
    }

    for (i=0; i<nelems; i++) buf[i] = 0;
    err = ncmpi_get_att_text(ncid, NC_GLOBAL, name, buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) {
        char expect = 'a' + i % 16;
        if (buf[i] != 'a' + i % 16) {
            printf("Error at %s line %d: expecting attr[%zd] value %c but got %c\n",
                   __FILE__,__LINE__,i,expect,buf[i]);
            nerrs++;
            break;
        }
    }

    name = "large_attr_1";
    err = ncmpi_inq_attlen(ncid, NC_GLOBAL, name, &inq_nelems);
    if (inq_nelems != nelems) {
        printf("Error at %s line %d: expecting attr %s nelems %lld but got %lld\n",
               __FILE__,__LINE__,name, nelems,inq_nelems);
        nerrs++;
    }

    for (i=0; i<nelems; i++) buf[i] = 0;
    err = ncmpi_get_att_text(ncid, NC_GLOBAL, name, buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) {
        char expect = 'a' + i % 16;
        if (buf[i] != 'a' + i % 16) {
            printf("Error at %s line %d: expecting attr[%zd] value %c but got %c\n",
                   __FILE__,__LINE__,i,expect,buf[i]);
            nerrs++;
            break;
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR
    if (nerrs > 0) goto err_out;

    /* create a new file and put 2 local attributes, total size > 2 GiB -----*/
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid);
    CHECK_ERR

    err = ncmpi_def_var(ncid, "var", NC_INT, 1, &dimid, &varid);
    CHECK_ERR

    for (i=0; i<nelems; i++) buf[i] = 'a' + i % 16;

    /* put two local attributes (total size > 2GiB) */
    name = "large_attr_0";
    err = ncmpi_put_att_text(ncid, varid, name, nelems, buf);
    if (!(cmode & NC_64BIT_DATA)) EXP_ERR(NC_EINVAL)
    else CHECK_ERR

    name = "large_attr_1";
    err = ncmpi_put_att_text(ncid, varid, name, nelems, buf);
    if (!(cmode & NC_64BIT_DATA)) EXP_ERR(NC_EINVAL)
    else CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    if (nerrs > 0) goto err_out;

    /* open the file and read back the two local attributes -----------------*/
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, info, &ncid);
    CHECK_ERR
    if (err != NC_NOERR) goto err_out;

    err = ncmpi_inq_varid(ncid, "var", &varid);
    CHECK_ERR

    name = "large_attr_0";
    err = ncmpi_inq_attlen(ncid, varid, name, &inq_nelems);
    if (inq_nelems != nelems) {
        printf("Error at %s line %d: expecting attr %s len %lld but got %lld\n",
               __FILE__,__LINE__,name,nelems,inq_nelems);
        nerrs++;
    }

    for (i=0; i<nelems; i++) buf[i] = 0;
    err = ncmpi_get_att_text(ncid, varid, name, buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) {
        char expect = 'a' + i % 16;
        if (buf[i] != 'a' + i % 16) {
            printf("Error at %s line %d: expecting attr[%zd] value %c but got %c\n",
                   __FILE__,__LINE__,i,expect,buf[i]);
            nerrs++;
            break;
        }
    }

    name = "large_attr_1";
    err = ncmpi_inq_attlen(ncid, varid, name, &inq_nelems);
    if (inq_nelems != nelems) {
        printf("Error at %s line %d: expecting attr %s len %lld but got %lld\n",
               __FILE__,__LINE__,name,nelems,inq_nelems);
        nerrs++;
    }

    for (i=0; i<nelems; i++) buf[i] = 0;
    err = ncmpi_get_att_text(ncid, varid, name, buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) {
        char expect = 'a' + i % 16;
        if (buf[i] != 'a' + i % 16) {
            printf("Error at %s line %d: expecting attr[%zd] value %c but got %c\n",
                   __FILE__,__LINE__,i,expect,buf[i]);
            nerrs++;
            break;
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR
    if (nerrs > 0) goto err_out;

    free(buf);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0) {
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
            ncmpi_inq_malloc_list();
        }
    }

    if (verbose) {
        err = ncmpi_inq_malloc_max_size(&malloc_size);
        printf("\n%d: PnetCDF internal memory footprint high water mark %.2f MB\n",
               rank, (float)malloc_size/1048576);
    }

err_out:
    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

