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
 * https://docs.unidata.ucar.edu/nug/current/md_types.html#data_type
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
int test_cdf12(const char *filename, int format, int coll_io, MPI_Info info)
{
    char val[MPI_MAX_INFO_VAL];
    int err, nerrs=0, ncid, vid, fvid, dvid, dimid, flag, bb_enabled=0;
    unsigned char uc;
    signed char sc;
    int si;
    double dbuf;

    /* check whether burst buffering is enabled */
    MPI_Info_get(info, "nc_burst_buf", MPI_MAX_INFO_VAL - 1, val, &flag);
    if (flag && strcasecmp(val, "enable") == 0) {
        bb_enabled = 1;
        /* does not work for NetCDF4 files when burst-buffering is enabled */
        if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC)
            return 0;
    }

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid); CHECK_ERR

    /* for CDF-1 and CDF-2, a special case is made: there is no NC_ERANGE
     * error can occur converting between NC_BYTE and unsigned char.
     * https://docs.unidata.ucar.edu/nug/current/md_types.html#data_type
     * In brief, NC_BYTE is signed in all signed CDF-2 APIs, and unsigned in
     * all unsigned APIs. In CDF-2, there is only one unsigned API, _uchar.
     */
    uc = 255;
    err = ncmpi_put_att_uchar(ncid, NC_GLOBAL, "att1", NC_BYTE, 1, &uc); CHECK_ERR
    uc = 0; /* initialize with a number that is not 0 */
    err = ncmpi_get_att_uchar(ncid, NC_GLOBAL, "att1", &uc); CHECK_ERR
    if (uc != 255) {
        printf("Error at line %d: unexpected read value %d (expecting 255)\n",__LINE__,(int)uc);
        assert(0);
    }
    sc = 3; /* initialize with a number that is not -1 or -0 */
    /* No NC_ERANGE as the internal and external types are considered the same */
    err = ncmpi_get_att_schar(ncid, NC_GLOBAL, "att1", &sc); CHECK_ERR
    if (   sc != -1     /* 2-complement bit representation */
        && sc != -0) {  /* 1-complement bit representation */
        printf("Error at line %d: unexpected read value %d (expecting 255)\n",__LINE__,(int)sc);
        assert(0);
    }

    /* expect NC_ERANGE */
    dbuf=NC_MAX_DOUBLE/2.0;
    err = ncmpi_put_att_double(ncid, NC_GLOBAL, "attf", NC_FLOAT, 1, &dbuf);
    EXP_ERR(NC_ERANGE)

    /* add a new attribute of NC_DOUBLE with a value > NC_MAX_FLOAT */
    err = ncmpi_put_att_double(ncid, NC_GLOBAL, "attd", NC_DOUBLE, 1, &dbuf);
    CHECK_ERR

#if defined(PNETCDF_ERANGE_FILL) && PNETCDF_ERANGE_FILL == 1
    if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
        float fbuf;
        /* read back attf and expect NC_FILL_FLOAT */
        dbuf = 0.0;
        err = ncmpi_get_att_double(ncid, NC_GLOBAL, "attf", &dbuf); CHECK_ERR
        if (dbuf != NC_FILL_FLOAT) {
            printf("Error at line %d: unexpected read value %f (expecting %f)\n",__LINE__,dbuf,NC_FILL_FLOAT);
            assert(0);
        }
        /* read back attd as float and expect NC_ERANGE */
        err = ncmpi_get_att_float(ncid, NC_GLOBAL, "attd", &fbuf); EXP_ERR(NC_ERANGE)
    }
#endif

    err = ncmpi_def_dim(ncid, "x", 1, &dimid); CHECK_ERR

    /* NC_BYTE is an 8-bit signed integer, i.e. C type of signed char */
    err = ncmpi_def_var(ncid, "var_byte", NC_BYTE, 1, &dimid, &vid);
    CHECK_ERR
    err = ncmpi_def_var(ncid, "var_flt", NC_FLOAT, 1, &dimid, &fvid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_dbl", NC_DOUBLE, 1, &dimid, &dvid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    /* No NC_ERANGE should be returned for CDF-1 and 2 */
    /* range of unsigned char is from 0 to 255. PnetCDF should detect the value
     * being out of range. But because NetCDF standard makes a special case for
     * CDF-1 and CDF-2, which treats external NC_BYTE as unsigned char with no
     * type conversion.
     *
     * https://docs.unidata.ucar.edu/nug/current/md_types.html#data_type
     * The _uchar and _schar functions were introduced in netCDF-3 to eliminate
     * an ambiguity, and support both signed and unsigned byte data. In
     * netCDF-2, whether the external NC_BYTE type represented signed or
     * unsigned values was left up to the user. In netcdf-3, we treat NC_BYTE
     * as signed for the purposes of conversion to short, int, long, float, or
     * double. (Of course, no conversion takes place when the internal type is
     * signed char.) In the _uchar functions, we treat NC_BYTE as if it were
     * unsigned. Thus, no NC_ERANGE error can occur converting between NC_BYTE
     * and unsigned char.
     */
    uc = 255;
    if (coll_io)
        err = ncmpi_put_var_uchar_all(ncid, vid, &uc);
    else
        err = ncmpi_put_var_uchar(ncid, vid, &uc);
    CHECK_ERR

    uc = 3; /* set read buffer to a different value before read */
    if (coll_io)
        err = ncmpi_get_var_uchar_all(ncid, vid, &uc);
    else
        err = ncmpi_get_var_uchar(ncid, vid, &uc);
    CHECK_ERR

    /* In netCDF-2, whether the external NC_BYTE type represented signed or
     * unsigned values was left up to the user.
     */
    if (format != NC_FORMAT_CLASSIC && format != NC_FORMAT_64BIT_OFFSET && uc != 255) {
        printf("Error %s at %d: uc[0] expect 255 but got %d\n",
                basename(__FILE__),__LINE__,(int)uc);
        assert(0);
    }

    /* No NC_ERANGE should be returned for CDF-1 and 2 */
    sc = -128;
    if (coll_io)
        err = ncmpi_put_var_schar_all(ncid, vid, &sc);
    else
        err = ncmpi_put_var_schar(ncid, vid, &sc);
    CHECK_ERR

    sc = 0; /* set read buffer to a different value before read */
    if (coll_io)
        err = ncmpi_get_var_schar_all(ncid, vid, &sc);
    else
        err = ncmpi_get_var_schar(ncid, vid, &sc);
    CHECK_ERR

    if (format != NC_FORMAT_CLASSIC && format != NC_FORMAT_64BIT_OFFSET && sc != -128) {
        printf("Error %s at %d: sc expect -128 but got %d\n",
                basename(__FILE__),__LINE__,(int)sc);
        assert(0);
    }

    /* expect NC_ERANGE */
    si = -129;
    if (coll_io)
        err = ncmpi_put_var_int_all(ncid, vid, &si);
    else
        err = ncmpi_put_var_int(ncid, vid, &si);
    if (bb_enabled) {
        CHECK_ERR
        err = ncmpi_flush(ncid);
    }
    EXP_ERR(NC_ERANGE)

    if (si != -129) { /* check if put buffer content is altered */
        printf("Error %s at %d: si expect -129 but got %d\n",
                basename(__FILE__),__LINE__,si);
        assert(0);
    }

    /* expect NC_ERANGE */
    si = 256;
    if (coll_io)
        err = ncmpi_put_var_int_all(ncid, vid, &si);
    else
        err = ncmpi_put_var_int(ncid, vid, &si);
    if (bb_enabled) {
        CHECK_ERR
        err = ncmpi_flush(ncid);
    }
    EXP_ERR(NC_ERANGE)

    if (si != 256) { /* check if put buffer content is altered */
        printf("Error %s at %d: si expect 256 but got %d\n",
                basename(__FILE__),__LINE__,si);
        assert(0);
    }

    /* expect no error */
    si = -128;
    if (coll_io)
        err = ncmpi_put_var_int_all(ncid, vid, &si);
    else
        err = ncmpi_put_var_int(ncid, vid, &si);
    CHECK_ERR

    si = 0; /* set read buffer to a different value before read */
    if (coll_io)
        err = ncmpi_get_var_int_all(ncid, vid, &si);
    else
        err = ncmpi_get_var_int(ncid, vid, &si);
    CHECK_ERR

    if (format != NC_FORMAT_CLASSIC && format != NC_FORMAT_64BIT_OFFSET && si != -128) {
        printf("Error %s at %d: si expect -128 but got %d\n",
                basename(__FILE__),__LINE__,si);
        assert(0);
    }

    /* expect NC_ERANGE */
    dbuf = NC_MAX_DOUBLE/2.0;
    if (coll_io)
        err = ncmpi_put_var_double_all(ncid, fvid, &dbuf);
    else
        err = ncmpi_put_var_double(ncid, fvid, &dbuf);
    if (bb_enabled) {
        CHECK_ERR
        err = ncmpi_flush(ncid);
    }
    EXP_ERR(NC_ERANGE)

    /* write a value > NC_MAX_FLOAT */
    if (coll_io)
        err = ncmpi_put_var_double_all(ncid, dvid, &dbuf);
    else
        err = ncmpi_put_var_double(ncid, dvid, &dbuf);
    CHECK_ERR


#if defined(PNETCDF_ERANGE_FILL) && PNETCDF_ERANGE_FILL == 1
    if (format != NC_FORMAT_NETCDF4 && format != NC_FORMAT_NETCDF4_CLASSIC) {
        float fbuf;
        /* read back and expect NC_FILL_FLOAT */
        dbuf = 0.0; /* set read buffer to a different value before read */
        if (coll_io)
            err = ncmpi_get_var_double_all(ncid, fvid, &dbuf);
        else
            err = ncmpi_get_var_double(ncid, fvid, &dbuf);
        CHECK_ERR
        if (dbuf != NC_FILL_FLOAT) {
            printf("Error %s at %d: dbuf expect NC_FILL_FLOAT (%f) but got %f\n",
                    basename(__FILE__),__LINE__,NC_FILL_FLOAT,dbuf);
            assert(0);
        }

        /* read back dvid as float and expect NC_ERANGE */
        fbuf = 0.0; /* set read buffer to a different value before read */
        if (coll_io)
            err = ncmpi_get_var_float_all(ncid, dvid, &fbuf);
        else
            err = ncmpi_get_var_float(ncid, dvid, &fbuf);
        EXP_ERR(NC_ERANGE)
    }
#endif

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

static
int test_cdf345(const char *filename, int format, int coll_io, MPI_Info info)
{
    char val[MPI_MAX_INFO_VAL];
    int err, nerrs=0, ncid, uc_vid, sc_vid, dimid, flag, bb_enabled=0;
    unsigned char uc;
    signed char sc;

    /* check whether burst buffering is enabled */
    MPI_Info_get(info, "nc_burst_buf", MPI_MAX_INFO_VAL - 1, val, &flag);
    if (flag && strcasecmp(val, "enable") == 0) {
        bb_enabled = 1;
        /* does not work for NetCDF4 files when burst-buffering is enabled */
        if (format == NC_FORMAT_NETCDF4 || format == NC_FORMAT_NETCDF4_CLASSIC)
            return 0;
    }

    /* Set format. */
    err = ncmpi_set_default_format(format, NULL);
    CHECK_ERR

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, info, &ncid); CHECK_ERR

    /* CDF-5 considers NC_BYTE a signed 1-byte integer and NC_UBYTE an
     * unsigned 1-byte integer. The special case in CDF-2 for skipping
     * NC_ERANGE checking for converting between NC_BYTE and unsigned
     * char is no longer held.
     */
    uc = 255;
    err = ncmpi_put_att_uchar(ncid, NC_GLOBAL, "att1", NC_UBYTE, 1, &uc); CHECK_ERR

    /* in CDF-5, get 255 to a schar buffer should result in NC_ERANGE */
    err = ncmpi_get_att_schar(ncid, NC_GLOBAL, "att1", &sc); EXP_ERR(NC_ERANGE)

    sc = -1; /* a value should cause NC_ERANGE */
    err = ncmpi_put_att_schar(ncid, NC_GLOBAL, "att2", NC_UBYTE, 1, &sc); EXP_ERR(NC_ERANGE)

    err = ncmpi_def_dim(ncid, "x", 1, &dimid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_ubyte", NC_UBYTE, 1, &dimid, &uc_vid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_byte",  NC_BYTE,  1, &dimid, &sc_vid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (!coll_io) {
        err = ncmpi_begin_indep_data(ncid);
        CHECK_ERR
    }

    uc = 255;
    if (coll_io)
        err = ncmpi_put_var_uchar_all(ncid, uc_vid, &uc);
    else
        err = ncmpi_put_var_uchar(ncid, uc_vid, &uc);
    CHECK_ERR

    /* in CDF-5, get 255 to an schar should result in NC_ERANGE */
    if (coll_io)
        err = ncmpi_get_var_schar_all(ncid, uc_vid, &sc);
    else
        err = ncmpi_get_var_schar(ncid, uc_vid, &sc);
    EXP_ERR(NC_ERANGE)

    sc = -1; /* in CDF-5, put -1 to an uchar should result in NC_ERANGE */
    if (coll_io)
        err = ncmpi_put_var_schar_all(ncid, uc_vid, &sc);
    else
        err = ncmpi_put_var_schar(ncid, uc_vid, &sc);
    if (bb_enabled) {
        CHECK_ERR
        err = ncmpi_flush(ncid);
    }
    EXP_ERR(NC_ERANGE)

    uc = 255; /* in CDF-5, put 255 to a schar should result in NC_ERANGE */
    if (coll_io)
        err = ncmpi_put_var_uchar_all(ncid, sc_vid, &uc);
    else
        err = ncmpi_put_var_uchar(ncid, sc_vid, &uc);
    if (bb_enabled) {
        CHECK_ERR
        err = ncmpi_flush(ncid);
    }
    EXP_ERR(NC_ERANGE)

    sc = -1;
    if (coll_io)
        err = ncmpi_put_var_schar_all(ncid, sc_vid, &sc);
    else
        err = ncmpi_put_var_schar(ncid, sc_vid, &sc);
    CHECK_ERR

    uc = 0; /* in CDF-5, get -1 to an uchar should result in NC_ERANGE */
    if (coll_io)
        err = ncmpi_get_var_uchar_all(ncid, sc_vid, &uc);
    else
        err = ncmpi_get_var_uchar(ncid, sc_vid, &uc);
    EXP_ERR(NC_ERANGE)

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

static
int test_io(const char *out_path,
            const char *in_path, /* ignored */
            int         format,
            int         coll_io,
            MPI_Info    info)
{
    int nerrs=0;

    if (format == NC_FORMAT_CLASSIC ||
        format == NC_FORMAT_64BIT_OFFSET ||
        format == NC_FORMAT_NETCDF4_CLASSIC)
        nerrs += test_cdf12(out_path, format, coll_io, info);
    else
        nerrs += test_cdf345(out_path, format, coll_io, info);

    return nerrs;
}

int main(int argc, char **argv) {

    int err;
    loop_opts opt;

    MPI_Init(&argc, &argv);

    opt.num_fmts = sizeof(nc_formats) / sizeof(int);
    opt.formats  = nc_formats;
    opt.ina      = 1; /* test intra-node aggregation */
    opt.drv      = 1; /* test PNCIO driver */
    opt.ind      = 1; /* test hint romio_no_indep_rw */
    opt.chk      = 0; /* test hint pnc_data_move_chunk_size */
    opt.bb       = 1; /* test burst-buffering feature */
    opt.mod      = 1; /* test independent data mode */
    opt.hdr_diff = 1; /* run ncmpidiff for file header only */
    opt.var_diff = 1; /* run ncmpidiff for variables */

    err = tst_main(argc, argv, "NC_ERANGE", opt, test_io);

    MPI_Finalize();

    return err;
}
