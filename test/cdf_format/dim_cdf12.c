/*********************************************************************
 *
 *  Copyright (C) 2016, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests CDF-1, CDF-2 file formats using the allowable maximal
 * dimension size. It also tests the following large-file support limitations.
 *
 * For CDF-1
 *
 * If there is no record variables defined in the file, only one fixed-size
 * variable can exceed 2 GiB in size (it can be as large as the underlying file
 * system permits.) It must be the last variable defined in the file, and the
 * offset to the beginning of this variable must be less than about 2 GiB.
 *
 * The limit is really 2^31 - 4. If you were to specify a variable size of 2^31
 * -3, for example, it would be rounded up to the nearest multiple of 4 bytes,
 * which would be 2^31, which is larger than the largest signed integer, 2^31
 * - 1.
 *
 * If you use the unlimited dimension, record variables may exceed 2 GiB in
 * size, as long as the offset of the start of each record variable within a
 * record is less than 2 GiB - 4.
 *
 *
 *
 * For CDF-2
 *
 * No fixed-size variable can require more than 2^32 - 4 bytes (i.e. 4GiB - 4
 * bytes, or 4,294,967,292 bytes) of storage for its data, unless it is the
 * last fixed-size variable and there are no record variables. When there are
 * no record variables, the last fixed-size variable can be any size supported
 * by the file system, e.g. terabytes.
 *
 * A 64-bit offset format netCDF file can have up to 2^32 - 1 fixed-size
 * variables, each under 4GiB in size. If there are no record variables in the
 * file the last fixed-size variable can be any size.
 *
 * No record variable can require more than 2^32 - 4 bytes of storage for each
 * record's worth of data, unless it is the last record variable. A 64-bit
 * offset format netCDF file can have up to 2^32 - 1 records, of up to 2^32 - 1
 * variables, as long as the size of one record's data for each record variable
 * except the last is less than 4 GiB - 4.
 *
 * Note also that all netCDF variables and records are padded to 4 byte boundaries.
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
            int         format,  /* ignored */
            int         coll_io, /* ignored */
            MPI_Info    info)
{
    int rank, nprocs, err, nerrs=0;
    int ncid, cmode, varid, dimid[3];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Note this test program must use the 512-byte alignment setting */

    /* use the 512-byte fixed-size variable starting file offset alignment,
     * which is also the header extent align size
     */
    MPI_Info_set(info, "nc_var_align_size", "512");

    /* create a new CDF-1 file ----------------------------------------------*/
    cmode = NC_CLOBBER;

    /* max dimension size for CDF-1 file is NC_MAX_INT */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", (MPI_Offset)1+NC_MAX_INT, &dimid[0]);
    EXP_ERR(NC_EDIMSIZE)
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT, &dimid[0]); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* use the max dimension size to define a 1D variable */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT, &dimid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_CHAR, 1, dimid, &varid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* use the max dimension size to define a 1D variable, followed by
     * another variable to make the file size > NC_MAX_INT */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 2,          &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_CHAR, 1, &dimid[0], &varid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT,  1, &dimid[1], &varid); CHECK_ERR
    /* for cdf-1, adding var1 after var0 will cause NC_EVARSIZE */
    err = ncmpi_close(ncid);
    EXP_ERR(NC_EVARSIZE)

    /* use the max dimension size - 1024 to define a 1D variable, followed
     * by another variable to make the file size < 2147483647 */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT-1024, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 2,               &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_CHAR, 1, &dimid[0], &varid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT,  1, &dimid[1], &varid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* define the first variable of type short that makes the file size >
     * NC_MAX_INT. error should be reported in ncmpi_enddef() or
     * ncmpi_close() */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 2,          &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_SHORT, 1, &dimid[0], &varid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_CHAR,  1, &dimid[1], &varid); CHECK_ERR
    err = ncmpi_close(ncid);
    EXP_ERR(NC_EVARSIZE)

    /* define two variables to make the file size just < 2147483647 */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT-512-8, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 2,                &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_CHAR, 1, &dimid[0], &varid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT,  1, &dimid[1], &varid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* define two variables to make the file size just > NC_MAX_INT */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT/2+1, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 2,              &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_INT, 1, &dimid[0], &varid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, &dimid[1], &varid); CHECK_ERR
    err = ncmpi_close(ncid);
    EXP_ERR(NC_EVARSIZE)

    /* create a new CDF-2 file ----------------------------------------------*/
    cmode = NC_CLOBBER | NC_64BIT_OFFSET;

    /* max dimension size for CDF-2 file is NC_MAX_INT */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", (MPI_Offset)1+NC_MAX_INT, &dimid[0]);
    EXP_ERR(NC_EDIMSIZE)
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT, &dimid[0]); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* use the max dimension size to define a 1D variable */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT, &dimid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_CHAR, 1, dimid, &varid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* use the max dimension size to define a 1D variable, followed by
     * another variable to make the file size > 2 * NC_MAX_INT */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 2,          &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_CHAR, 1, &dimid[0], &varid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT,  1, &dimid[1], &varid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* define the first variable of type short that makes the file size >
     * 4294967295. error should be reported in ncmpi_enddef() or
     * ncmpi_close() */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 2,       &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_SHORT, 1, &dimid[0], &varid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT,   1, &dimid[1], &varid); CHECK_ERR
    err = ncmpi_close(ncid);
    EXP_ERR(NC_EVARSIZE)

    /* define 2 1D int variables of dimension size > NC_MAX_INT */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 2,          &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_INT, 1, &dimid[0], &varid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, &dimid[1], &varid); CHECK_ERR
    err = ncmpi_close(ncid);
    EXP_ERR(NC_EVARSIZE)

    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT/2+1, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 2,              &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_INT, 1, &dimid[0], &varid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, &dimid[1], &varid); CHECK_ERR
    err = ncmpi_close(ncid);
    EXP_ERR(NC_EVARSIZE)

    /* No record variable can require more than 2^32 - 4 bytes of storage for
     * each record's worth of data, unless it is the last record variable.
     */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Z", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT/64,   &dimid[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 64,              &dimid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_INT, 3, dimid, &varid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 3, dimid, &varid); CHECK_ERR
    err = ncmpi_close(ncid);
    EXP_ERR(NC_EVARSIZE)

    /* test large record variable that is not defined last */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Z", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT/64,   &dimid[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 64,              &dimid[2]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_INT, 3, dimid, &varid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 2, dimid, &varid); CHECK_ERR
    err = ncmpi_close(ncid);
    EXP_ERR(NC_EVARSIZE)

    /* Note for developers: keep the last test that produces no error, so the
     * output file can be tested by ncvalidator in wrap_runs.sh
     */
    err = ncmpi_create(MPI_COMM_WORLD, out_path, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", NC_MAX_INT/2, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 2,            &dimid[1]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var0", NC_INT, 1, &dimid[0], &varid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, &dimid[1], &varid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR
    err = ncmpi_open(MPI_COMM_WORLD, out_path, NC_NOWRITE, info, &ncid); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

/*----< main() >--------------------------------------------------------------*/
int main(int argc, char **argv) {

    int err, formats[] = {0};
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(formats) / sizeof(int);
    opt.formats  = formats;
    opt.ina      = 0; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 0; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 0; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "defining dim in CDF-1/2 format", opt, test_io);

    MPI_Finalize();

    return err;
}

