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
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

static
int test_cdf2(char *filename)
{
    int err, nerrs=0, ncid, vid, dimid;
    unsigned char uc[1];
    signed char sc[1];
    int si[1];
#ifdef BUILD_DRIVER_DW
    int dw_enabled=0;
#endif

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid); CHECK_ERR

#ifdef BUILD_DRIVER_DW
    {
        int flag;
        char hint[MPI_MAX_INFO_VAL];
        MPI_Info infoused;
        ncmpi_inq_file_info(ncid, &infoused);
        MPI_Info_get(infoused, "nc_dw", MPI_MAX_INFO_VAL - 1, hint, &flag);
        if (flag && strcasecmp(hint, "enable") == 0)
            dw_enabled = 1;
        MPI_Info_free(&infoused);
    }
#endif

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

    err = ncmpi_def_dim(ncid, "x", 1, &dimid); CHECK_ERR
    err = ncmpi_def_var(ncid, "var_byte", NC_BYTE, 1, &dimid, &vid); CHECK_ERR
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
#ifdef BUILD_DRIVER_DW
    if (dw_enabled) {
        CHECK_ERR
        err = ncmpi_sync(ncid);
    }
#endif
    EXP_ERR(NC_ERANGE)

    if (si[0] != -129) { /* check if put buffer content is altered */
        printf("Error at line %d: put buffer content altered %d (expecting -128)\n",__LINE__,si[0]);
        nerrs++;
    }

    /* expect NC_ERANGE */
    si[0] = 256;
    err = ncmpi_put_var_int_all(ncid, vid, si);
#ifdef BUILD_DRIVER_DW
    if (dw_enabled) {
        CHECK_ERR
        err = ncmpi_sync(ncid);
    }
#endif
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

    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

static
int test_cdf5(char *filename)
{
    int err, nerrs=0, ncid, uc_vid, sc_vid, dimid;
    unsigned char uc[1];
    signed char sc[1];
#ifdef BUILD_DRIVER_DW
    int dw_enabled=0;
#endif

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER|NC_64BIT_DATA, MPI_INFO_NULL, &ncid); CHECK_ERR

#ifdef BUILD_DRIVER_DW
    {
        int flag;
        char hint[MPI_MAX_INFO_VAL];
        MPI_Info infoused;
        ncmpi_inq_file_info(ncid, &infoused);
        MPI_Info_get(infoused, "nc_dw", MPI_MAX_INFO_VAL - 1, hint, &flag);
        if (flag && strcasecmp(hint, "enable") == 0)
            dw_enabled = 1;
        MPI_Info_free(&infoused);
    }
#endif

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
#ifdef BUILD_DRIVER_DW
    if (dw_enabled){
        CHECK_ERR
        err = ncmpi_sync(ncid);
    }
#endif
    EXP_ERR(NC_ERANGE)

    uc[0] = 255; /* in CDF-5, put 255 to a schar should result in NC_ERANGE */
    err = ncmpi_put_var_uchar_all(ncid, sc_vid, uc);
#ifdef BUILD_DRIVER_DW
    if (dw_enabled){
        CHECK_ERR
        err = ncmpi_sync(ncid);
    }
#endif
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
    char filename[256];
    int err, nerrs=0, rank;

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

    nerrs += test_cdf2(filename);
    nerrs += test_cdf5(filename);

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}
