/*
 *  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id$
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests conflicted in-memory data types against external data
 * types and see if the correct error codes were returned.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o check_type check_type.c -lpnetcdf
 *
 *    % mpiexec -l -n 1 check_type testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <pnetcdf.h>

#include <testutils.h>

#define EXP_ERR_N_TYPE(expect_err,itype,etype) { \
    if (err != expect_err) { \
        printf("Error at line %d in %s: itype=%9s etype=%-9s err=%s\n", \
               __LINE__,__FILE__,itype,etype_name(etype),ncmpi_strerrno(err)); \
        nerrs++; \
    }  \
}

static char* etype_name(nc_type etype) {
    switch (etype) {
        case (0):  return "NC_NAT";
        case (1):  return "NC_BYTE";
        case (2):  return "NC_CHAR";
        case (3):  return "NC_SHORT";
        case (4):  return "NC_INT";
        case (5):  return "NC_FLOAT";
        case (6):  return "NC_DOUBLE";
        case (7):  return "NC_UBYTE";
        case (8):  return "NC_USHORT";
        case (9):  return "NC_UINT";
        case (10): return "NC_INT64";
        case (11): return "NC_UINT64";
        default:
              return "Invalid nc_type";
    }
}

int main(int argc, char* argv[])
{
    char filename[256], *varname[12], buf[1024], attname[256];
    int i, err, nerrs=0, rank, ncid, dimid, varid[12];

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
        sprintf(cmd_str, "*** TESTING C   %s for checking for type conflict ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }

    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER|NC_64BIT_DATA, MPI_INFO_NULL, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "x", 2, &dimid); CHECK_ERR

    varname[0]  = "var_nat";
    varname[1]  = "var_byte";
    varname[2]  = "var_char";
    varname[3]  = "var_short";
    varname[4]  = "var_int";
    varname[5]  = "var_float";
    varname[6]  = "var_double";
    varname[7]  = "var_ubyte";
    varname[8]  = "var_ushort";
    varname[9]  = "var_uint";
    varname[10] = "var_int64";
    varname[11] = "var_uint64";
    for (i=NC_BYTE; i<=NC_UINT64; i++) {
        err = ncmpi_def_var(ncid, varname[i], i, 1, &dimid, &varid[i]); CHECK_ERR
    }

    for (i=0; i<1024; i++) buf[i]=i;
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "att_text", 3, (char*)buf); CHECK_ERR
    EXP_ERR_N_TYPE(NC_NOERR, "text", NC_CHAR)

    for (i=0; i<1024; i++) buf[i]=0;
    for (i=NC_BYTE; i<=NC_UINT64; i++) {
        int expect_err = NC_NOERR;

        if (i == NC_CHAR) expect_err = NC_ECHAR;

        sprintf(attname,"att_uchar_for_var_%s",varname[i]);
        err = ncmpi_put_att_uchar    (ncid, NC_GLOBAL, attname, i, 3, (unsigned char*)      buf);
        EXP_ERR_N_TYPE(expect_err, "uchar", i)
        sprintf(attname,"att_schar_for_var_%s",varname[i]);
        err = ncmpi_put_att_schar    (ncid, NC_GLOBAL, attname, i, 3, (signed char*)        buf);
        EXP_ERR_N_TYPE(expect_err, "schar", i)
        sprintf(attname,"att_short_for_var_%s",varname[i]);
        err = ncmpi_put_att_short    (ncid, NC_GLOBAL, attname, i, 3, (short*)              buf);
        EXP_ERR_N_TYPE(expect_err, "short", i)
        sprintf(attname,"att_int_for_var_%s",varname[i]);
        err = ncmpi_put_att_int      (ncid, NC_GLOBAL, attname, i, 3, (int*)                buf);
        EXP_ERR_N_TYPE(expect_err, "int", i)
        sprintf(attname,"att_float_for_var_%s",varname[i]);
        err = ncmpi_put_att_float    (ncid, NC_GLOBAL, attname, i, 3, (float*)              buf);
        EXP_ERR_N_TYPE(expect_err, "float", i)
        sprintf(attname,"att_double_for_var_%s",varname[i]);
        err = ncmpi_put_att_double   (ncid, NC_GLOBAL, attname, i, 3, (double*)             buf);
        EXP_ERR_N_TYPE(expect_err, "double", i)
        sprintf(attname,"att_ushort_for_var_%s",varname[i]);
        err = ncmpi_put_att_ushort   (ncid, NC_GLOBAL, attname, i, 3, (unsigned short*)     buf);
        EXP_ERR_N_TYPE(expect_err, "ushort", i)
        sprintf(attname,"att_uint_for_var_%s",varname[i]);
        err = ncmpi_put_att_uint     (ncid, NC_GLOBAL, attname, i, 3, (unsigned int*)       buf);
        EXP_ERR_N_TYPE(expect_err, "uint", i)
        sprintf(attname,"att_longlong_for_var_%s",varname[i]);
        err = ncmpi_put_att_longlong (ncid, NC_GLOBAL, attname, i, 3, (long long*)          buf);
        EXP_ERR_N_TYPE(expect_err, "longlong", i)
        sprintf(attname,"att_ulonglong_for_var_%s",varname[i]);
        err = ncmpi_put_att_ulonglong(ncid, NC_GLOBAL, attname, i, 3, (unsigned long long*) buf);
        EXP_ERR_N_TYPE(expect_err, "ulonglong", i)
    }

    for (i=NC_BYTE; i<=NC_UINT64; i++) {
        int expect_err = NC_NOERR;

        if (i == NC_CHAR) continue;

        sprintf(attname,"att_uchar_for_var_%s",varname[i]);
        err = ncmpi_get_att_uchar    (ncid, NC_GLOBAL, attname, (unsigned char*)      buf);
        EXP_ERR_N_TYPE(expect_err, "uchar", i)
        sprintf(attname,"att_schar_for_var_%s",varname[i]);
        err = ncmpi_get_att_schar    (ncid, NC_GLOBAL, attname, (signed char*)        buf);
        EXP_ERR_N_TYPE(expect_err, "schar", i)
        sprintf(attname,"att_short_for_var_%s",varname[i]);
        err = ncmpi_get_att_short    (ncid, NC_GLOBAL, attname, (short*)              buf);
        EXP_ERR_N_TYPE(expect_err, "short", i)
        sprintf(attname,"att_int_for_var_%s",varname[i]);
        err = ncmpi_get_att_int      (ncid, NC_GLOBAL, attname, (int*)                buf);
        EXP_ERR_N_TYPE(expect_err, "int", i)
        sprintf(attname,"att_float_for_var_%s",varname[i]);
        err = ncmpi_get_att_float    (ncid, NC_GLOBAL, attname, (float*)              buf);
        EXP_ERR_N_TYPE(expect_err, "float", i)
        sprintf(attname,"att_double_for_var_%s",varname[i]);
        err = ncmpi_get_att_double   (ncid, NC_GLOBAL, attname, (double*)             buf);
        EXP_ERR_N_TYPE(expect_err, "double", i)
        sprintf(attname,"att_ushort_for_var_%s",varname[i]);
        err = ncmpi_get_att_ushort   (ncid, NC_GLOBAL, attname, (unsigned short*)     buf);
        EXP_ERR_N_TYPE(expect_err, "ushort", i)
        sprintf(attname,"att_uint_for_var_%s",varname[i]);
        err = ncmpi_get_att_uint     (ncid, NC_GLOBAL, attname, (unsigned int*)       buf);
        EXP_ERR_N_TYPE(expect_err, "uint", i)
        sprintf(attname,"att_longlong_for_var_%s",varname[i]);
        err = ncmpi_get_att_longlong (ncid, NC_GLOBAL, attname, (long long*)          buf);
        EXP_ERR_N_TYPE(expect_err, "longlong", i)
        sprintf(attname,"att_ulonglong_for_var_%s",varname[i]);
        err = ncmpi_get_att_ulonglong(ncid, NC_GLOBAL, attname, (unsigned long long*) buf);
        EXP_ERR_N_TYPE(expect_err, "ulonglong", i)
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
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();

    return (nerrs > 0);
}

