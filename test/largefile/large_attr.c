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

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    char *name, *buf;
    size_t i;
    int err, nerrs=0, ncid, varid, dimid;
    MPI_Offset nelems, inq_nelems;
    int rank, nprocs, color;
    MPI_Comm comm=MPI_COMM_NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    color = 1;

    if (nprocs > 2) {
        /* run on 2 ranks only, as this test allocates memory > 4GB per rank */
        /* split MPI_COMM_WORLD based on 'color' and use the same rank order */
        color = (rank < 2) ? 1 : 0;
        MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);
    }
    else
        comm = MPI_COMM_WORLD;

    if (!color) goto err_out;

    nelems = (MPI_Offset)NC_MAX_INT + 17;
    buf = (char*) malloc(nelems);

    /* Set file format */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    /* create a new file and put a large global attribute -------------------*/
    err = ncmpi_create(comm, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* put a large (> 2GiB) global attribute */
    for (i=0; i<nelems; i++) buf[i] = 'a' + i % 16;
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "large_attr", nelems, buf);
    if (format != NC_FORMAT_64BIT_DATA) EXP_ERR(NC_EINVAL)
    else CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* open the file and read back the large global attribute ---------------*/
    err = ncmpi_open(comm, out_path, NC_NOWRITE, info, &ncid);
    CHECK_ERR

    err = ncmpi_inq_attlen(ncid, NC_GLOBAL, "large_attr", &inq_nelems);
    CHECK_ERR
    if (inq_nelems != nelems) {
        printf("Error at %s line %d: expecting attr nelems "OFFFMT" but got "OFFFMT"\n",
               __FILE__,__LINE__,nelems,inq_nelems);
        nerrs++;
        goto err_out;
    }

    for (i=0; i<nelems; i++) buf[i] = 0;
    err = ncmpi_get_att_text(ncid, NC_GLOBAL, "large_attr", buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) {
        char expect = 'a' + i % 16;
        if (buf[i] != expect) {
            printf("Error at %s line %d: expecting attr[%zd] value %c but got %c\n",
                   __FILE__,__LINE__,i,expect,buf[i]);
            nerrs++;
            goto err_out;
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR
    if (nerrs > 0) goto err_out;

    /* create a new file and put a large local attribute -------------------*/
    err = ncmpi_create(comm, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid);
    CHECK_ERR

    err = ncmpi_def_var(ncid, "var", NC_INT, 1, &dimid, &varid);
    CHECK_ERR

    /* put a large (> 2GiB) global attribute */
    for (i=0; i<nelems; i++) buf[i] = 'a' + i % 16;
    err = ncmpi_put_att_text(ncid, varid, "large_attr", nelems, buf);
    if (format != NC_FORMAT_64BIT_DATA) EXP_ERR(NC_EINVAL)
    else CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* open the file and read back the large global attribute ---------------*/
    err = ncmpi_open(comm, out_path, NC_NOWRITE, info, &ncid);
    CHECK_ERR
    if (err != NC_NOERR) goto err_out;

    err = ncmpi_inq_varid(ncid, "var", &varid);
    CHECK_ERR

    err = ncmpi_inq_attlen(ncid, varid, "large_attr", &inq_nelems);
    if (inq_nelems != nelems) {
        printf("Error at %s line %d: expecting attr len "OFFFMT" but got "OFFFMT"\n",
               __FILE__,__LINE__,nelems,inq_nelems);
        nerrs++;
        goto err_out;
    }

    for (i=0; i<nelems; i++) buf[i] = 0;
    err = ncmpi_get_att_text(ncid, varid, "large_attr", buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) {
        char expect = 'a' + i % 16;
        if (buf[i] != expect) {
            printf("Error at %s line %d: expecting attr[%zd] value %c but got %c\n",
                   __FILE__,__LINE__,i,expect,buf[i]);
            nerrs++;
            goto err_out;
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR

    /* create a new file and put 2 global attributes, total size > 2 GiB ----*/
    nelems /= 2;

    err = ncmpi_create(comm, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    /* put two global attributes (total size > 2GiB) */
    for (i=0; i<nelems; i++) buf[i] = 'a' + i % 16;
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "large_attr_0", nelems, buf);
    if (format != NC_FORMAT_64BIT_DATA) EXP_ERR(NC_EINVAL)
    else CHECK_ERR

    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "large_attr_1", nelems, buf);
    if (format != NC_FORMAT_64BIT_DATA) EXP_ERR(NC_EINVAL)
    else CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* open the file and read back the large global attributes --------------*/
    err = ncmpi_open(comm, out_path, NC_NOWRITE, info, &ncid);
    CHECK_ERR
    if (err != NC_NOERR) goto err_out;

    name = "large_attr_0";
    err = ncmpi_inq_attlen(ncid, NC_GLOBAL, name, &inq_nelems);
    if (inq_nelems != nelems) {
        printf("Error at %s line %d: expecting attr %s nelems "OFFFMT" but got "OFFFMT"\n",
               __FILE__,__LINE__,name, nelems,inq_nelems);
        nerrs++;
        goto err_out;
    }

    for (i=0; i<nelems; i++) buf[i] = 0;
    err = ncmpi_get_att_text(ncid, NC_GLOBAL, name, buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) {
        char expect = 'a' + i % 16;
        if (buf[i] != expect) {
            printf("Error at %s line %d: expecting attr[%zd] value %c but got %c\n",
                   __FILE__,__LINE__,i,expect,buf[i]);
            nerrs++;
            goto err_out;
        }
    }

    name = "large_attr_1";
    err = ncmpi_inq_attlen(ncid, NC_GLOBAL, name, &inq_nelems);
    if (inq_nelems != nelems) {
        printf("Error at %s line %d: expecting attr %s nelems "OFFFMT" but got "OFFFMT"\n",
               __FILE__,__LINE__,name, nelems,inq_nelems);
        nerrs++;
        goto err_out;
    }

    for (i=0; i<nelems; i++) buf[i] = 0;
    err = ncmpi_get_att_text(ncid, NC_GLOBAL, name, buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) {
        char expect = 'a' + i % 16;
        if (buf[i] != expect) {
            printf("Error at %s line %d: expecting attr[%zd] value %c but got %c\n",
                   __FILE__,__LINE__,i,expect,buf[i]);
            nerrs++;
            goto err_out;
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR

    /* create a new file and put 2 local attributes, total size > 2 GiB -----*/
    err = ncmpi_create(comm, out_path, NC_CLOBBER, info, &ncid);
    CHECK_ERR

    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid);
    CHECK_ERR

    err = ncmpi_def_var(ncid, "var", NC_INT, 1, &dimid, &varid);
    CHECK_ERR

    for (i=0; i<nelems; i++) buf[i] = 'a' + i % 16;

    /* put two local attributes (total size > 2GiB) */
    name = "large_attr_0";
    err = ncmpi_put_att_text(ncid, varid, name, nelems, buf);
    if (format != NC_FORMAT_64BIT_DATA) EXP_ERR(NC_EINVAL)
    else CHECK_ERR

    name = "large_attr_1";
    err = ncmpi_put_att_text(ncid, varid, name, nelems, buf);
    if (format != NC_FORMAT_64BIT_DATA) EXP_ERR(NC_EINVAL)
    else CHECK_ERR

    err = ncmpi_enddef(ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* open the file and read back the two local attributes -----------------*/
    err = ncmpi_open(comm, out_path, NC_NOWRITE, info, &ncid);
    CHECK_ERR
    if (err != NC_NOERR) goto err_out;

    err = ncmpi_inq_varid(ncid, "var", &varid);
    CHECK_ERR

    name = "large_attr_0";
    err = ncmpi_inq_attlen(ncid, varid, name, &inq_nelems);
    if (inq_nelems != nelems) {
        printf("Error at %s line %d: expecting attr %s len "OFFFMT" but got "OFFFMT"\n",
               __FILE__,__LINE__,name,nelems,inq_nelems);
        nerrs++;
        goto err_out;
    }

    for (i=0; i<nelems; i++) buf[i] = 0;
    err = ncmpi_get_att_text(ncid, varid, name, buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) {
        char expect = 'a' + i % 16;
        if (buf[i] != expect) {
            printf("Error at %s line %d: expecting attr[%zd] value %c but got %c\n",
                   __FILE__,__LINE__,i,expect,buf[i]);
            nerrs++;
            goto err_out;
        }
    }

    name = "large_attr_1";
    err = ncmpi_inq_attlen(ncid, varid, name, &inq_nelems);
    if (inq_nelems != nelems) {
        printf("Error at %s line %d: expecting attr %s len "OFFFMT" but got "OFFFMT"\n",
               __FILE__,__LINE__,name,nelems,inq_nelems);
        nerrs++;
        goto err_out;
    }

    for (i=0; i<nelems; i++) buf[i] = 0;
    err = ncmpi_get_att_text(ncid, varid, name, buf);
    CHECK_ERR

    for (i=0; i<nelems; i++) {
        char expect = 'a' + i % 16;
        if (buf[i] != expect) {
            printf("Error at %s line %d: expecting attr[%zd] value %c but got %c\n",
                   __FILE__,__LINE__,i,expect,buf[i]);
            nerrs++;
            goto err_out;
        }
    }

    err = ncmpi_close(ncid); CHECK_ERR

    free(buf);

err_out:
    if (comm != MPI_COMM_WORLD && comm != MPI_COMM_NULL)
        MPI_Comm_free(&comm);

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    int formats[] = {NC_FORMAT_64BIT_DATA};
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 0; /* test intra-node aggregation */
    opt.drv      = 0; /* test PNCIO driver */
    opt.ind      = 0; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint nc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 0; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "one large ATTR", opt, test_io);

    MPI_Finalize();

    return err;
}
