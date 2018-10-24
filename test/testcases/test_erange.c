/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests whether NC_ERANGE error code can be reported correctly.
 * Note in CDF-1 and CDF-2, a special case is made to NOT report NC_ERANGE
 * when the variable is of NC_BYTE type and the calling APIs are of uchar. See
 * http://www.unidata.ucar.edu/software/netcdf/docs/data_type.html#type_conversion
 *
 * In CDF-5, NC_ERANGE is checked for when the external data type mismatches the
 * internal one.
 *
 * The test uses the following 2 case.
 * 1. get a value of 255 from a NC_UBYTE variable defined in a netCDF file to a
 *    memory buffer of signed char through e.g. API ncmpi_get_var_schar_all
 * 2. put a value of -1 of signed char from an in-memory buffer to a NC_UBYTE
 *    variable defined in a netCDF file
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o test_erange test_erange.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 test_erange testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

static
int test_cdf12(char *filename, int bb_enabled, int cmode)
{
    int err, nerrs=0, ncid, vid, fvid, dvid, dimid;
    unsigned char uc[1];
    signed char sc[1];
    int si[1];
    double dbuf;

    cmode |= NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid); CHECK_ERR

    /* for CDF-1 and CDF-2, a special case is made: there is no NC_ERANGE
     * error can occur converting between NC_BYTE and unsigned char.
     * http://www.unidata.ucar.edu/software/netcdf/docs/data_type.html#type_conversion
     * In brief, NC_BYTE is signed in all signed CDF-2 APIs, and unsigned in
     * all unsigned APIs. In CDF-2, there is only one unsigned API, _uchar.
     */
    uc[0] = 255;
    err = ncmpi_put_att_uchar(ncid, NC_GLOBAL, "att1", NC_BYTE, 1, uc); CHECK_ERR
    uc[0] = 0; /* initialize with a number that is not 0 */
    err = ncmpi_get_att_uchar(ncid, NC_GLOBAL, "att1", uc); CHECK_ERR
    if (uc[0] != 255) {
        printf("Error at line %d: unexpected read value %d (expecting 255)\n",__LINE__,(int)uc[0]);
        nerrs++;
    }
    sc[0] = 3; /* initialize with a number that is not -1 or -0 */
    /* No NC_ERANGE as the internal and external types are considered the same */
    err = ncmpi_get_att_schar(ncid, NC_GLOBAL, "att1", sc); CHECK_ERR
    if (   sc[0] != -1     /* 2-complement bit representation */
        && sc[0] != -0) {  /* 1-complement bit representation */
        printf("Error at line %d: unexpected read value %d (expecting 255)\n",__LINE__,(int)uc[0]);
        nerrs++;
    }

    /* expect NC_ERANGE */
    dbuf=NC_MAX_DOUBLE/2.0;
    err = ncmpi_put_att_double(ncid, NC_GLOBAL, "attf", NC_FLOAT, 1, &dbuf);
    EXP_ERR(NC_ERANGE)

    /* add a new attribute of NC_DOUBLE with a value > NC_MAX_FLOAT */
    err = ncmpi_put_att_double(ncid, NC_GLOBAL, "attd", NC_DOUBLE, 1, &dbuf);
    CHECK_ERR

#if defined(PNETCDF_ERANGE_FILL) && PNETCDF_ERANGE_FILL == 1
    if (! (cmode & NC_NETCDF4)) {
        float fbuf;
        /* read back attf and expect NC_FILL_FLOAT */
        dbuf = 0.0;
        err = ncmpi_get_att_double(ncid, NC_GLOBAL, "attf", &dbuf); CHECK_ERR
        if (dbuf != NC_FILL_FLOAT) {
            printf("Error at line %d: unexpected read value %f (expecting %f)\n",__LINE__,dbuf,NC_FILL_FLOAT);
            nerrs++;
        }
        /* read back attd as float and expect NC_ERANGE */
        err = ncmpi_get_att_float(ncid, NC_GLOBAL, "attd", &fbuf); EXP_ERR(NC_ERANGE)
    }
#endif

    err = ncmpi_def_dim(ncid, "x", 1, &dimid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_byte", NC_BYTE, 1, &dimid, &vid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_flt", NC_FLOAT, 1, &dimid, &fvid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_dbl", NC_DOUBLE, 1, &dimid, &dvid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* No NC_ERANGE should be returned for CDF-1 and 2 */
    uc[0] = 255;
    err = ncmpi_put_var_uchar_all(ncid, vid, uc); CHECK_ERR
    uc[0] = 3; /* initialize with a number that is not -1 or -0 */
    err = ncmpi_get_var_uchar_all(ncid, vid, uc); CHECK_ERR
    if (uc[0] != 255) {
        printf("Error at line %d: unexpected read value %d (expecting 255)\n",__LINE__,(int)uc[0]);
        nerrs++;
    }

    /* No NC_ERANGE should be returned for CDF-1 and 2 */
    sc[0] = -128;
    err = ncmpi_put_var_schar_all(ncid, vid, sc); CHECK_ERR
    sc[0] = 0;
    err = ncmpi_get_var_schar_all(ncid, vid, sc); CHECK_ERR
    if (sc[0] != -128) {
        printf("Error at line %d: unexpected read value %d (expecting -128)\n",__LINE__,(int)sc[0]);
        nerrs++;
    }

    /* expect NC_ERANGE */
    si[0] = -129;
    err = ncmpi_put_var_int_all(ncid, vid, si);
    if (bb_enabled) {
        CHECK_ERR
        err = ncmpi_flush(ncid);
    }
    EXP_ERR(NC_ERANGE)

    if (si[0] != -129) { /* check if put buffer content is altered */
        printf("Error at line %d: put buffer content altered %d (expecting -128)\n",__LINE__,si[0]);
        nerrs++;
    }

    /* expect NC_ERANGE */
    si[0] = 256;
    err = ncmpi_put_var_int_all(ncid, vid, si);
    if (bb_enabled) {
        CHECK_ERR
        err = ncmpi_flush(ncid);
    }
    EXP_ERR(NC_ERANGE)

    if (si[0] != 256) { /* check if put buffer content is altered */
        printf("Error at line %d: put buffer content altered %d (expecting 256)\n",__LINE__,si[0]);
        nerrs++;
    }

    /* expect no error */
    si[0] = -128;
    err = ncmpi_put_var_int_all(ncid, vid, si); CHECK_ERR
    si[0] = 0;
    err = ncmpi_get_var_int_all(ncid, vid, si); CHECK_ERR
    if (si[0] != -128) {
        printf("Error at line %d: unexpected read value %d (expecting -128)\n",__LINE__,si[0]);
        nerrs++;
    }

    /* expect NC_ERANGE */
    dbuf = NC_MAX_DOUBLE/2.0;
    err = ncmpi_put_var_double_all(ncid, fvid, &dbuf);
    if (bb_enabled) {
        CHECK_ERR
        err = ncmpi_flush(ncid);
    }
    EXP_ERR(NC_ERANGE)

    /* write a value > NC_MAX_FLOAT */
    err = ncmpi_put_var_double_all(ncid, dvid, &dbuf); CHECK_ERR

#if defined(PNETCDF_ERANGE_FILL) && PNETCDF_ERANGE_FILL == 1
    if (! (cmode & NC_NETCDF4)) {
        float fbuf;
        /* read back and expect NC_FILL_FLOAT */
        dbuf = 0.0;
        err = ncmpi_get_var_double_all(ncid, fvid, &dbuf); CHECK_ERR
        if (dbuf != NC_FILL_FLOAT) {
            printf("Error at line %d: unexpected read value %f (expecting %f)\n",__LINE__,dbuf,NC_FILL_FLOAT);
            nerrs++;
        }

        /* read back dvid as float and expect NC_ERANGE */
        err = ncmpi_get_var_float_all(ncid, dvid, &fbuf); EXP_ERR(NC_ERANGE)
    }
#endif

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

static
int test_cdf345(char *filename, int bb_enabled, int cmode)
{
    int err, nerrs=0, ncid, uc_vid, sc_vid, dimid;
    unsigned char uc[1];
    signed char sc[1];

    cmode |= NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid); CHECK_ERR

    /* CDF-5 considers NC_BYTE a signed 1-byte integer and NC_UBYTE an
     * unsigned 1-byte integer. The special case in CDF-2 for skipping
     * NC_ERANGE checking for converting between NC_BYTE and unsigned
     * char is no longer held.
     */
    uc[0] = 255;
    err = ncmpi_put_att_uchar(ncid, NC_GLOBAL, "att1", NC_UBYTE, 1, uc); CHECK_ERR

    /* in CDF-5, get 255 to a schar buffer should result in NC_ERANGE */
    err = ncmpi_get_att_schar(ncid, NC_GLOBAL, "att1", sc); EXP_ERR(NC_ERANGE)

    sc[0] = -1; /* a value should cause NC_ERANGE */
    err = ncmpi_put_att_schar(ncid, NC_GLOBAL, "att2", NC_UBYTE, 1, sc); EXP_ERR(NC_ERANGE)

    err = ncmpi_def_dim(ncid, "x", 1, &dimid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_ubyte", NC_UBYTE, 1, &dimid, &uc_vid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_byte",  NC_BYTE,  1, &dimid, &sc_vid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    uc[0] = 255;
    err = ncmpi_put_var_uchar_all(ncid, uc_vid, uc); CHECK_ERR

    /* in CDF-5, get 255 to an schar should result in NC_ERANGE */
    err = ncmpi_get_var_schar_all(ncid, uc_vid, sc); EXP_ERR(NC_ERANGE)

    sc[0] = -1; /* in CDF-5, put -1 to an uchar should result in NC_ERANGE */
    err = ncmpi_put_var_schar_all(ncid, uc_vid, sc);
    if (bb_enabled) {
        CHECK_ERR
        err = ncmpi_flush(ncid);
    }
    EXP_ERR(NC_ERANGE)

    uc[0] = 255; /* in CDF-5, put 255 to a schar should result in NC_ERANGE */
    err = ncmpi_put_var_uchar_all(ncid, sc_vid, uc);
    if (bb_enabled) {
        CHECK_ERR
        err = ncmpi_flush(ncid);
    }
    EXP_ERR(NC_ERANGE)

    sc[0] = -1;
    err = ncmpi_put_var_schar_all(ncid, sc_vid, sc); CHECK_ERR
    uc[0] = 0; /* in CDF-5, get -1 to an uchar should result in NC_ERANGE */
    err = ncmpi_get_var_uchar_all(ncid, sc_vid, uc); EXP_ERR(NC_ERANGE)

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

int main(int argc, char* argv[])
{
    char filename[256], *hint_value;
    int err, nerrs=0, rank, bb_enabled=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
        sprintf(cmd_str, "*** TESTING C   %s for checking for NC_ERANGE ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    /* check whether burst buffering is enabled */
    if (inq_env_hint("nc_burst_buf", &hint_value)) {
        if (strcasecmp(hint_value, "enable") == 0) bb_enabled = 1;
        free(hint_value);
    }

    nerrs += test_cdf12(filename, bb_enabled, 0);
    nerrs += test_cdf12(filename, bb_enabled, NC_64BIT_OFFSET);
#if ENABLE_NETCDF4
    if (!bb_enabled) {
        nerrs += test_cdf12(filename, bb_enabled, NC_NETCDF4 | NC_CLASSIC_MODEL);
        nerrs += test_cdf345(filename, bb_enabled, NC_NETCDF4);
    }
#endif
    nerrs += test_cdf345(filename, bb_enabled, NC_64BIT_DATA);

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
