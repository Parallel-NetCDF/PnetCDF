/*********************************************************************
 *
 *  Copyright (C) 2012, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

/*    This example is similar to collective_write.c but using nonblocking APIs,
 *    iput or bput (default is iput). It creates a netcdf file in CD-5 format
 *    and writes a number of 3D integer non-record variables. The measured
 *    write bandwidth is reported at the end.
 *
 *    To compile:
 *        mpicc -O2 nonblocking_write.c -o nonblocking_write -lpnetcdf
 *
 *    To run:
 *        mpiexec -n num_processes ./nonblocking_write -l len file_name
 *    where len decides the size of each local array, which is len x len x len.
 *    Each global variable is of size len*len*len * nprocs * sizeof(int).
 *    All variables are partitioned among all processes in a 3D
 *    block-block-block fashion. Below is an example standard output from
 *    command:
 *        mpiexec -n 32 ./nonblocking_write /pvfs2/wkliao/testfile.nc -l 100
 *
 *    MPI hint: cb_nodes        = 2
 *    MPI hint: cb_buffer_size  = 16777216
 *    MPI hint: striping_factor = 32
 *    MPI hint: striping_unit   = 1048576
 *    Local array size 100 x 100 x 100 integers, size = 3.81 MB
 *    Global array size 400 x 400 x 200 integers, write size = 0.30 GB
 *     procs    Global array size  exec(sec)  write(MB/s)
 *     -------  ------------------  ---------  -----------
 *        32     400 x  400 x  200     6.67       45.72
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <mpi.h>
#include <pnetcdf.h>

#define NTIMES     2
#define NDIMS      3
#define SCA_NVARS  5
#define FIX_NVARS  5
#define REC_NVARS 10

static int verbose;

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h | -q | -c | -b | -l len] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [-c] Allocate all write buffers in a contiguous space (default no)\n"
    "       [-b] Use bput APIs instead of iput APIs (default iput)\n"
    "       [-l len] size of each dimension of the local array\n"
    "       [file_name] output netCDF file name\n";
    fprintf(stderr, help, argv0);
}

/*----< pnetcdf_check_mem_usage() >------------------------------------------*/
/* check PnetCDF library internal memory usage */
static int
pnetcdf_check_mem_usage(MPI_Comm comm)
{
    int err, nerrs=0, rank;
    MPI_Offset malloc_size, sum_size;

    MPI_Comm_rank(comm, &rank);

    /* print info about PnetCDF internal malloc usage */
    err = ncmpi_inq_malloc_max_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && verbose)
            printf("maximum heap memory allocated by PnetCDF internally is %lld bytes\n",
                   sum_size);

        /* check if there is any PnetCDF internal malloc residue */
        err = ncmpi_inq_malloc_size(&malloc_size);
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }
    else if (err != NC_ENOTENABLED) {
        printf("Error at %s:%d: %s\n", __FILE__,__LINE__,ncmpi_strerror(err));
        nerrs++;
    }
    return nerrs;
}

/*----< print_info() >------------------------------------------------------*/
static
void print_info(MPI_Info *info_used)
{
    int     flag;
    char    info_cb_nodes[64], info_cb_buffer_size[64];
    char    info_striping_factor[64], info_striping_unit[64];

    strcpy(info_cb_nodes,        "undefined");
    strcpy(info_cb_buffer_size,  "undefined");
    strcpy(info_striping_factor, "undefined");
    strcpy(info_striping_unit,   "undefined");

    MPI_Info_get(*info_used, "cb_nodes", 64, info_cb_nodes, &flag);
    MPI_Info_get(*info_used, "cb_buffer_size", 64, info_cb_buffer_size, &flag);
    MPI_Info_get(*info_used, "striping_factor", 64, info_striping_factor, &flag);
    MPI_Info_get(*info_used, "striping_unit", 64, info_striping_unit, &flag);

    printf("MPI hint: cb_nodes        = %s\n", info_cb_nodes);
    printf("MPI hint: cb_buffer_size  = %s\n", info_cb_buffer_size);
    printf("MPI hint: striping_factor = %s\n", info_striping_factor);
    printf("MPI hint: striping_unit   = %s\n", info_striping_unit);
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    extern int optind;
    extern char *optarg;
    int i, j, k, err, nerrs=0, debug=0, use_contig_buf=0, use_bput=0;
    int nprocs, len=0, nelems, rank;
    int *sca_buf, *fix_buf[FIX_NVARS], *rec_buf[REC_NVARS];
    int gsizes[NDIMS], psizes[NDIMS];
    double write_timing, max_write_timing, write_bw;
    char filename[256], str[512];
    int ncid, dimids[NDIMS+1];
    int sca_var[SCA_NVARS], fix_var[FIX_NVARS], rec_var[REC_NVARS], *req, *st;
    MPI_Offset starts[NDIMS+1], counts[NDIMS+1], write_size, sum_write_size;
    MPI_Offset header_size, bbufsize, put_size, sum_put_size;
    MPI_Info info, info_used;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    verbose = 1;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hqcbl:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'c': use_contig_buf = 1;
                      break;
            case 'b': use_bput = 1;
                      break;
            case 'l': len = atoi(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

    len = (len <= 0) ? 10 : len;

    for (i=0; i<NDIMS; i++) psizes[i] = 0;

    MPI_Dims_create(nprocs, NDIMS, psizes);
    if (rank == 0 && verbose && debug)
        printf("Process DIMS psizes=%2d %2d %2d\n",
               psizes[0],psizes[1],psizes[2]);

    starts[0] = 0;
    starts[1] = (rank / (psizes[1] * psizes[2])) % psizes[0];
    starts[2] = (rank / psizes[2]) % psizes[1];
    starts[3] = rank % psizes[2];

    counts[0] = 1;
    nelems = 1;
    for (i=0; i<NDIMS; i++) {
        gsizes[i]    = len * psizes[i];
        nelems      *= len;
        starts[i+1] *= len;
        counts[i+1]  = len;
    }

    if (verbose && debug)
        printf("%2d: starts=%2lld %2lld %2lld %2lld counts=%2lld %2lld %2lld %2lld\n",
               rank, starts[0], starts[1], starts[2], starts[3],
                     counts[0], counts[1], counts[2], counts[3]);

    /* allocate buffers */
    if (use_contig_buf) {
        /* all write buffers are allocated in a single contiguous space */
        size_t total_nelems;

        total_nelems = SCA_NVARS + (FIX_NVARS + REC_NVARS) * nelems;

        sca_buf = (int*) malloc(total_nelems * sizeof(int));

        for (i=0; i<FIX_NVARS; i++)
            fix_buf[i] = sca_buf + SCA_NVARS + nelems * i;

        for (i=0; i<REC_NVARS; i++)
            rec_buf[i] = sca_buf + SCA_NVARS + FIX_NVARS * nelems + nelems * i;
    }
    else {
        /* allocate individual buffers separately +1 ensure non-contiguity*/
        sca_buf = (int*) malloc((SCA_NVARS+1) * sizeof(int));
        for (i=0; i<FIX_NVARS; i++)
            fix_buf[i] = (int *) malloc((nelems+1) * sizeof(int));
        for (i=0; i<REC_NVARS; i++)
            rec_buf[i] = (int *) malloc((nelems+1) * sizeof(int));
    }

    /* initialize buffer contents */
    for (j=0; j<SCA_NVARS; j++) sca_buf[j] = rank + j;
    for (i=0; i<FIX_NVARS; i++) {
        for (j=0; j<nelems; j++) fix_buf[i][j] = rank;
    }
    for (i=0; i<REC_NVARS; i++) {
        for (j=0; j<nelems; j++) rec_buf[i][j] = rank;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    write_timing = MPI_Wtime();

    /* set an MPI-IO hint to disable file offset alignment for fixed-size
     * variables */
    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_var_align_size", "1");

    /* disable PnetCDF internal buffering */
    MPI_Info_set(info, "nc_ibuf_size", "0");
    MPI_Info_set(info, "nc_in_place_swap", "enable");

    /* create the file */
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER|NC_64BIT_DATA,
                       info, &ncid);
    if (err != NC_NOERR) {
        printf("Error at line %d in %s: ncmpi_create() file %s (%s)\n",
        __LINE__,__FILE__,filename,ncmpi_strerror(err));
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(1);
    }

    MPI_Info_free(&info);

    req = (int*) malloc((SCA_NVARS + FIX_NVARS + REC_NVARS) * sizeof(int));
    st  = (int*) malloc((SCA_NVARS + FIX_NVARS + REC_NVARS) * sizeof(int));

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimids[0]);
    ERR
    for (i=0; i<NDIMS; i++) {
        sprintf(str, "%c", 'z'-i);
        err = ncmpi_def_dim(ncid, str, gsizes[i], &dimids[i+1]);
        ERR
    }

    /* define scalar variables */
    for (i=0; i<SCA_NVARS; i++) {
        sprintf(str, "scalar_var_%d", i);
        err = ncmpi_def_var(ncid, str, NC_INT, 0, NULL, &sca_var[i]);
        ERR
    }

    /* define fix-sized variables */
    for (i=0; i<FIX_NVARS; i++) {
        sprintf(str, "fix_var_%d", i);
        err = ncmpi_def_var(ncid, str, NC_INT, NDIMS, dimids+1, &fix_var[i]);
        ERR
    }

    /* define record variables */
    for (i=0; i<REC_NVARS; i++) {
        sprintf(str, "rec_var_%d", i);
        err = ncmpi_def_var(ncid, str, NC_INT, NDIMS+1, dimids, &rec_var[i]);
        ERR
    }

    /* exit the define mode */
    err = ncmpi_enddef(ncid);
    ERR

    /* get all the MPI-IO and PnetCDF hints used */
    err = ncmpi_inq_file_info(ncid, &info_used);
    ERR

    if (!use_bput) {
        k = 0;

        /* only rank 0 writes scalar variables using iput */
        if (rank == 0) {
            for (i=0; i<SCA_NVARS; i++) {
                err = ncmpi_iput_var_int(ncid, sca_var[i], &sca_buf[i], &req[k++]);
                ERR
            }
        }

        /* write one fix-sized variable at a time using iput */
        for (i=0; i<FIX_NVARS; i++) {
            err = ncmpi_iput_vara_int(ncid, fix_var[i], starts+1, counts+1, fix_buf[i], &req[k++]);
            ERR
        }

        /* write one record variable at a time using iput */
        for (j=0; j<NTIMES; j++) {
            starts[0] = j;
            for (i=0; i<REC_NVARS; i++) {
                err = ncmpi_iput_vara_int(ncid, rec_var[i], starts, counts, rec_buf[i], &req[k++]);
                ERR
            }

            /* wait for the nonblocking I/O to complete */
            err = ncmpi_wait_all(ncid, k, req, st);
            ERR
            for (i=0; i<k; i++) {
                if (st[i] != NC_NOERR)
                    printf("Error at line %d in %s: nonblocking write fails on request %d (%s)\n",
                    __LINE__,__FILE__,i, ncmpi_strerror(st[i]));
            }
            k = 0;
        }
    }
    else { /* write using bput */

        /* bbufsize must be max of data type converted before and after */
        bbufsize = (SCA_NVARS + nelems * (FIX_NVARS + REC_NVARS)) * sizeof(int);
        err = ncmpi_buffer_attach(ncid, bbufsize);
        ERR

        k = 0;

        /* only rank 0 writes scalar variables time using bput */
        if (rank == 0) {
            for (i=0; i<SCA_NVARS; i++) {
                err = ncmpi_bput_var_int(ncid, sca_var[i], &sca_buf[i], &req[k++]);
                ERR
            }
        }

        /* write one fix_sized variable at a time using bput */
        for (i=0; i<FIX_NVARS; i++) {
            err = ncmpi_bput_vara_int(ncid, fix_var[i], starts+1, counts+1, fix_buf[i], &req[k++]);
            ERR
            /* can safely change contents of fix_buf[i] here */
        }

        /* write one record variable at a time using bput */
        for (j=0; j<NTIMES; j++) {
            starts[0] = j;
            for (i=0; i<REC_NVARS; i++) {
                err = ncmpi_bput_vara_int(ncid, rec_var[i], starts, counts, rec_buf[i], &req[k++]);
                ERR
            }

            /* wait for the nonblocking I/O to complete */
            err = ncmpi_wait_all(ncid, k, req, st);
            ERR
            for (i=0; i<k; i++) {
                if (st[i] != NC_NOERR)
                    printf("Error at line %d in %s: nonblocking write fails on request %d (%s)\n",
                    __LINE__,__FILE__,i, ncmpi_strerror(st[i]));
            }
            k = 0;
        }

        /* detach the temporary buffer */
        err = ncmpi_buffer_detach(ncid);
        ERR
    }

    ncmpi_inq_put_size(ncid, &put_size);
    ncmpi_inq_header_size(ncid, &header_size);

    /* close the file */
    err = ncmpi_close(ncid);
    ERR

    if (use_contig_buf)
        free(sca_buf);
    else {
        free(sca_buf);
        for (i=0; i<FIX_NVARS; i++) free(fix_buf[i]);
        for (i=0; i<REC_NVARS; i++) free(rec_buf[i]);
    }
    free(req);
    free(st);

    write_timing = MPI_Wtime() - write_timing;

    write_size = nelems * (FIX_NVARS + NTIMES * REC_NVARS) + SCA_NVARS;
    write_size *= sizeof(int);

    MPI_Reduce(&write_size,   &sum_write_size,   1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&put_size,     &sum_put_size,     1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&write_timing, &max_write_timing, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    /* rank 0 also writes header and updates the record number to the file
     * header NTIME times
     */
    if (rank == 0 && verbose) {
        sum_write_size += header_size + NTIMES * 8;

        printf("\n");
        if (use_bput)
            printf("Using PnetCDF nonblocking APIs: bput\n");
        else
            printf("Using PnetCDF nonblocking APIs: iput\n");
        if (use_contig_buf)
            printf("All write buffers are in a contiguous space\n");
        else
            printf("All write buffers are allocated separately\n");
        printf("Total amount writes                     (include header) = %lld bytes\n", sum_write_size);
        printf("Total amount writes reported by pnetcdf (include header) = %lld bytes\n", sum_put_size);
        printf("\n");
        float subarray_size = (float)nelems*sizeof(int)/1048576.0;
        print_info(&info_used);
        printf("Local array size %d x %d x %d integers, size = %.2f MB\n",len,len,len,subarray_size);
        sum_write_size /= 1048576.0;
        printf("Global array size %d x %d x %d integers, write size = %.2f GB\n",
               gsizes[0], gsizes[1], gsizes[2], sum_write_size/1024.0);

        write_bw = sum_write_size/max_write_timing;
        printf(" procs    Global array size  exec(sec)  write(MB/s)\n");
        printf("-------  ------------------  ---------  -----------\n");
        printf(" %4d    %4d x %4d x %4d %8.2f  %10.2f\n", nprocs,
               gsizes[0], gsizes[1], gsizes[2], max_write_timing, write_bw);
    }
    MPI_Info_free(&info_used);

    nerrs += pnetcdf_check_mem_usage(MPI_COMM_WORLD);

    MPI_Finalize();
    return (nerrs > 0);
}

