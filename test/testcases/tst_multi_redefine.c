/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This program tests entering define modes multiple times, growing header
 * extension causing moving data section to plances of higher offsets, and
 * checking the contents of all variables defined.
 *
 * The compile and run commands are given below.
 *
 *    % mpicc -g -o tst_multi_redefine tst_multi_redefine.c -lpnetcdf
 *
 *    % mpiexec -l -n 4 ./tst_multi_redefine testfile.nc
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp() */
#include <libgen.h>  /* basename() */
#include <unistd.h>  /* getopt() */

#include <pnetcdf.h>

#include <testutils.h>

#define NROUNDS 2
#define NY 10
#define NX 10

static int verbose;

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef long long longlong;
typedef unsigned long long ulonglong;

#define END_DEF { \
    err = ncmpi_enddef(ncid); CHECK_ERROUT \
    err = ncmpi_inq_header_size(ncid, &new_hdr_size); CHECK_ERROUT \
    err = ncmpi_inq_header_extent(ncid, &new_hdr_ext); CHECK_ERROUT \
    if (verbose && rank == 0) { \
        printf("Add var %2d: header      size grows from %4lld to %4lld\n", \
               k, old_hdr_size, new_hdr_size); \
        if (new_hdr_ext > old_hdr_ext) \
            printf("Add var %2d: header extension grows from %4lld to %4lld\n", \
                k, old_hdr_ext, new_hdr_ext); \
    } \
    old_hdr_size = new_hdr_size; \
    old_hdr_ext  = new_hdr_ext; \
}

#define EXP_VAL(nn) ((j + nn + rank) % NC_MAX_BYTE)

#define DEF_VAR(xtype) { \
    err = ncmpi_redef(ncid); CHECK_ERROUT \
    sprintf(str, "fix_var_%d", k); \
    err = ncmpi_def_var(ncid, str, xtype, 2, dimids+1, &fix_varids[k]); \
    CHECK_ERROUT \
    sprintf(str, "attribute of fix-sized variable %d", k); \
    err = ncmpi_put_att_text(ncid, fix_varids[k], "attr", strlen(str), str); \
    CHECK_ERROUT \
    sprintf(str, "rec_var_%d", k); \
    err = ncmpi_def_var(ncid, str, xtype, 3, dimids, &rec_varids[k]); \
    CHECK_ERROUT \
    sprintf(str, "attribute of record variable %d", k); \
    err = ncmpi_put_att_text(ncid, rec_varids[k], "attr", strlen(str), str); \
    CHECK_ERROUT \
    END_DEF \
}

#define PUT_BUF_CHAR { \
    char *buf = (char*) malloc(sizeof(char) * nelems); \
    for (j=0; j<nelems; j++) buf[j] = EXP_VAL(k); \
    err = ncmpi_put_vara_text_all(ncid, fix_varids[k], start+1, count+1, buf); \
    CHECK_ERROUT \
    err = ncmpi_put_vara_text_all(ncid, rec_varids[k], start, count, buf); \
    CHECK_ERROUT \
    free(buf); \
}

#define PUT_BUF(itype) { \
    itype *buf = (itype*) malloc(sizeof(itype) * nelems); \
    for (j=0; j<nelems; j++) buf[j] = EXP_VAL(k); \
    err = ncmpi_put_vara_##itype##_all(ncid, fix_varids[k], start+1, count+1, buf); \
    CHECK_ERROUT \
    err = ncmpi_put_vara_##itype##_all(ncid, rec_varids[k], start, count, buf); \
    CHECK_ERROUT \
    free(buf); \
}

#define CHECK_BUF_CHAR { \
    char *buf = (char*) malloc(sizeof(char) * nelems); \
    err = ncmpi_get_vara_text_all(ncid, fix_varids[j], start+1, count+1, buf); \
    CHECK_ERROUT \
    for (m=0; m<nelems; m++) { \
        char exp = EXP_VAL(m); \
        if (buf[m] != exp) { \
            printf("Error: var %d expect buf[%d] = %d but got %d\n", \
                   j, m, (int)exp, (int)buf[m]); \
            nerrs++; \
            goto err_out; \
        } \
    } \
    memset(buf, 0, sizeof(char) * nelems); \
    err = ncmpi_get_vara_text_all(ncid, rec_varids[j], start, count, buf); \
    CHECK_ERROUT \
    for (m=0; m<nelems; m++) { \
        char exp = EXP_VAL(m); \
        if (buf[m] != exp) { \
            printf("Error: var %d expect buf[%d] = %d but got %d\n", \
                   j, m, (int)exp, (int)buf[m]); \
            nerrs++; \
            goto err_out; \
        } \
    } \
    free(buf); \
}

#define CHECK_BUF \
    for (j=0; j<k; j++) { \
        if (j % ntypes == 0) { \
            CHECK_BUF_CHAR \
            continue; \
        } \
        double *buf = (double*) malloc(sizeof(double) * nelems); \
        err = ncmpi_get_vara_double_all(ncid, fix_varids[j], start+1, count+1, buf); \
        CHECK_ERROUT \
        for (m=0; m<nelems; m++) { \
            double exp = EXP_VAL(m); \
            if (buf[m] != exp) { \
                printf("Error: var %d expect buf[%d] = %.f but got %.f\n", \
                       j, m, exp, buf[m]); \
                nerrs++; \
                goto err_out; \
            } \
        } \
        memset(buf, 0, sizeof(double) * nelems); \
        err = ncmpi_get_vara_double_all(ncid, rec_varids[j], start, count, buf); \
        CHECK_ERROUT \
        for (m=0; m<nelems; m++) { \
            double exp = EXP_VAL(m); \
            if (buf[m] != exp) { \
                printf("Error: var %d expect buf[%d] = %.f but got %.f\n", \
                       j, m, exp, buf[m]); \
                nerrs++; \
                goto err_out; \
            } \
        } \
        free(buf); \
    }

static int
tst_fmt(char *filename,
        int   cmode,
        int   nRounds,
        int   len_y,
        int   len_x)
{
    char str[64];
    int i, j, k, m, rank, nprocs, ncid, err, nerrs=0, dimids[3];
    int *fix_varids=NULL, *rec_varids=NULL, ntypes;
    MPI_Offset old_hdr_size, new_hdr_size;
    MPI_Offset old_hdr_ext, new_hdr_ext;
    MPI_Offset nelems, start[3], count[3];

    MPI_Info info=MPI_INFO_NULL;

    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    if (verbose && rank == 0) {
        char *fmt = "NC_CLASSIC";
        if (cmode == NC_64BIT_OFFSET)
            fmt = "NC_64BIT_OFFSET";
        else if (cmode == NC_64BIT_DATA)
            fmt = "NC_64BIT_DATA";
        printf("\n---- Testing file format %s ----\n", fmt);
    }

    /* CDF-1 and 2 have 6 external data types, CDF-5 has 11 */
    ntypes = (cmode & NC_64BIT_DATA) ? 11 : 6;

    fix_varids = (int*) malloc(sizeof(int) * nRounds * ntypes);
    rec_varids = (int*) malloc(sizeof(int) * nRounds * ntypes);

    start[0] = 0;
    start[1] = 0;
    start[2] = rank * len_x;
    count[0] = 1;
    count[1] = len_y;
    count[2] = len_x;
    nelems = count[1] * count[2];

    /* create a new file */
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERROUT

    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimids[0]); CHECK_ERROUT
    err = ncmpi_def_dim(ncid, "Y",    len_y,        &dimids[1]); CHECK_ERROUT
    err = ncmpi_def_dim(ncid, "X",    len_x*nprocs, &dimids[2]); CHECK_ERROUT

    /* Default v_align of 512 is pretty bad choice, as almost every iteration
     * grows the file extension.
     */
    // err = ncmpi__enddef(ncid, 0, 4096, 0, 0);
    err = ncmpi_enddef(ncid);
    CHECK_ERROUT

    err = ncmpi_inq_header_size(ncid, &old_hdr_size); CHECK_ERROUT
    err = ncmpi_inq_header_extent(ncid, &old_hdr_ext); CHECK_ERROUT
    if (verbose && rank == 0)
        printf("Newly created file with header size %lld, extension %lld\n",
               old_hdr_size, old_hdr_ext);

    k = 0;
    for (i=0; i<nRounds; i++) {
        /* add 2 new variables of type NC_CHAR */
        DEF_VAR(NC_CHAR)
        PUT_BUF_CHAR
        k++;
        CHECK_BUF

        /* add 2 new variables of type NC_BYTE */
        DEF_VAR(NC_BYTE)
        PUT_BUF(schar)
        k++;
        CHECK_BUF

        /* add 2 new variables of type NC_SHORT */
        DEF_VAR(NC_SHORT)
        PUT_BUF(short)
        k++;
        CHECK_BUF

        /* add 2 new variables of type NC_INT */
        DEF_VAR(NC_INT)
        PUT_BUF(int)
        k++;
        CHECK_BUF

        /* add 2 new variables of type NC_FLOAT */
        DEF_VAR(NC_FLOAT)
        PUT_BUF(float)
        k++;
        CHECK_BUF

        /* add 2 new variables of type NC_DOUBLE */
        DEF_VAR(NC_DOUBLE)
        PUT_BUF(double)
        k++;
        CHECK_BUF

        if (cmode & NC_64BIT_DATA) {
            /* add 2 new variables of type NC_UBYTE */
            DEF_VAR(NC_UBYTE)
            PUT_BUF(uchar)
            k++;
            CHECK_BUF

            /* add 2 new variables of type NC_USHORT */
            DEF_VAR(NC_USHORT)
            PUT_BUF(ushort)
            k++;
            CHECK_BUF

            /* add 2 new variables of type NC_UINT */
            DEF_VAR(NC_UINT)
            PUT_BUF(uint)
            k++;
            CHECK_BUF

            /* add 2 new variables of type NC_INT64 */
            DEF_VAR(NC_INT64)
            PUT_BUF(longlong)
            k++;
            CHECK_BUF

            /* add 2 new variables of type NC_UINT64 */
            DEF_VAR(NC_UINT64)
            PUT_BUF(ulonglong)
            k++;
            CHECK_BUF
        }
    }

    err = ncmpi_close(ncid); CHECK_ERROUT

err_out:
    if (fix_varids != NULL) free(fix_varids);
    if (rec_varids != NULL) free(rec_varids);
    if (info != MPI_INFO_NULL) MPI_Info_free(&info);

    return nerrs;
}

#define FILE_NAME "testfile.nc"

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [OPTIONS]...[filename]\n"
    "       [-h] Print help\n"
    "       [-v]: enable verbose output mode (default: no)\n"
    "       [-n num]: number of rounds entering define mode (default: %d)\n"
    "       [-y num]: Y dimension size of global array (default: %d)\n"
    "       [-x num]: X dimension size of local array (default: %d)\n"
    "       [filename]: output netCDF file name (default: %s)\n";
    fprintf(stderr, help, argv0, FILE_NAME, NROUNDS, NY, NX);
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    char filename[256];
    int i, fmt, rank, err, nerrs=0, cmode[3], nRounds, len_y, len_x;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 0;
    nRounds = NROUNDS;
    len_y = NY;
    len_x = NX;
    while ((i = getopt(argc, argv, "hvn:y:x:")) != EOF)
        switch(i) {
            case 'v': verbose = 1;
                      break;
            case 'n': nRounds = atoi(optarg);
                      break;
            case 'y': len_y = atoi(optarg);
                      break;
            case 'x': len_x = atoi(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, FILE_NAME);
    else                      snprintf(filename, 256, "%s", argv[optind]);

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for entering define mode ", basename(argv[0]));
        printf("%-66s ------ ", cmd_str); fflush(stdout);
        free(cmd_str);
    }
    cmode[0] = 0;
    cmode[1] = NC_64BIT_OFFSET;
    cmode[2] = NC_64BIT_DATA;

    for (fmt=0; fmt<3; fmt++) {
        nerrs += tst_fmt(filename, cmode[fmt], nRounds, len_y, len_x);
        if (nerrs > 0) goto err_out;
    }

    /* check if PnetCDF freed all internal malloc */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has "OFFFMT" bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
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

