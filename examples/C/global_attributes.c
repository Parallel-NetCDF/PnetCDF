/*********************************************************************
 *
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example creates a new file and add 2 global attributes, one is of text
 * type and the other 10 short integers. The file is closed and re-opened to
 * read back the attributes.
 *
 *    To compile:
 *        mpicc -O2 global_attributes.c -o global_attributes -lpnetcdf
 *
 * Example commands for MPI run and outputs from running ncmpidump on the
 * netCDF file produced by this example program:
 *
 *    % mpiexec -n 4 ./global_attributes /pvfs2/wkliao/testfile.nc
 *
 *    % ncmpidump /pvfs2/wkliao/testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-1
 *
 *    // global attributes:
 *                    :history = "Sun May  4 13:11:47 2014\n",
 *        "" ;
 *                    :digits = 0s, 1s, 2s, 3s, 4s, 5s, 6s, 7s, 8s, 9s ;
 *    }
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy(), strlen() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <mpi.h>
#include <pnetcdf.h>

#define ERR { \
    if(err!=NC_NOERR) { \
        printf("Error at line %d in %s: %s\n", __LINE__,__FILE__, ncmpi_strerror(err)); \
        nerrs++; \
        goto fn_exit; \
    } \
}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [filename] output netCDF file name\n";
    fprintf(stderr, help, argv0);
}

int main(int argc, char** argv)
{
    extern int optind;
    char filename[256];
    char str_att[128], att_name[NC_MAX_NAME];
    int i, rank, err, nerrs=0, verbose=1, ncid, cmode, omode, ngatts;
    short short_att[10], digit[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    time_t ltime;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hq")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid);
    ERR

    /* add a global attribute named "history": a time stamp at rank 0 */
    ltime = time(NULL); /* get the current calendar time */
    asctime_r(localtime(&ltime), str_att);

    /* make sure the time string are consistent among all processes */
    MPI_Bcast(str_att, strlen(str_att), MPI_CHAR, 0, MPI_COMM_WORLD);

    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "history", strlen(str_att),
                             &str_att[0]);
    ERR
    if (rank == 0 && verbose)
        printf("writing global attribute \"history\" of text %s\n",
               str_att);

    /* add another global attribute named "digits": an array of short type */
    err = ncmpi_put_att_short(ncid, NC_GLOBAL, "digits", NC_SHORT, 10,
                             &digit[0]);
    ERR
    if (rank == 0 && verbose)
        printf("writing global attribute \"digits\" of 10 short integers\n");

    /* close file */
    err = ncmpi_close(ncid);
    ERR

    /* open the newly created file for read only -----------------------------*/
    omode = NC_NOWRITE;
    err = ncmpi_open(MPI_COMM_WORLD, filename, omode, MPI_INFO_NULL, &ncid);
    ERR

    /* find the number of global attributes */
    err = ncmpi_inq_natts(ncid, &ngatts);
    ERR

    err = 0;
    if (ngatts != 2) {
        printf("Error at line %d in %s: expected number of global attributes is 2, but got %d\n",
               __LINE__,__FILE__,ngatts);
        err = -1;
    }
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (err < 0) { nerrs++; goto fn_exit; }

    /* find the name of first global attribute */
    err = ncmpi_inq_attname(ncid, NC_GLOBAL, 0, att_name);
    ERR

    err = 0;
    if (strncmp(att_name, "history", strlen("history"))) {
        printf("Error at line %d in %s: expected attribute name \"history\", but got %s\n",
               __LINE__,__FILE__,att_name);
        err = -1;
    }
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (err < 0) { nerrs++; goto fn_exit; }

    /* read attribute value */
    err = ncmpi_get_att_text(ncid, NC_GLOBAL, att_name, &str_att[0]);
    ERR

    /* find the name of second global attribute */
    err = ncmpi_inq_attname(ncid, NC_GLOBAL, 1, att_name);
    ERR

    err = 0;
    if (strncmp(att_name, "digits", strlen("digits"))) {
        printf("Error at line %d in %s: expected attribute name \"digits\", but got %s\n",
               __LINE__,__FILE__,att_name);
        err = -1;
    }
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (err < 0) { nerrs++; goto fn_exit; }

    /* read attribute value */
    err = ncmpi_get_att_short(ncid, NC_GLOBAL, att_name, &short_att[0]);
    ERR

    /* close file */
    err = ncmpi_close(ncid);
    ERR

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

fn_exit:
    MPI_Finalize();
    return (nerrs > 0);
}

