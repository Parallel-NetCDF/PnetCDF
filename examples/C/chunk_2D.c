/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * This example shows how to use the chunking and compression features of
 * PnetCDF to write a 3D record variable of integer data type in parallel. It
 * first defines netCDF variables, each of size NTIMES x NY x NX, where NTIMES,
 * NY, and NX are predefined constant.
 *
 * The data partitioning pattern is a checkerboard style, along both Y and X
 * dimensions. Each process writes a subarray per time record.
 *
 *    To compile:
 *        mpicc -O2 chunk_2D.c -o chunk_2D  \
 *                  -I/path/to/PnetCDF/include \
 *                  -I/path/to/ZLIB/include \
 *                  -I/path/to/SZ/include \
 *                  -L/path/to/PnetCDF/lib \
 *                  -L/path/to/ZLIB/lib \
 *                  -L/path/to/SZ/lib \
 *                  -lpnetcdf -lz -lm -ldl -lSZ -lzstd
 *
 * Example commands for MPI run and outputs from running ncmpidump on the
 * output netCDF file produced by this example program:
 *
 *    % mpiexec -n 4 ./chunk_2D testfile.nc
 *
 *    % ncmpidump testfile.nc
 *    netcdf testfile {
 *    // file format: CDF-5 (big variables)
 *    dimensions:
 *        time = UNLIMITED ; // (0 currently)
 *        Y = 10 ;
 *        X = 10 ;
 *        _datablock_dim_0 = 131484 ;
 *        _datablock_dim_1 = 412 ;
 *    variables:
 *        int var_0 ;
 *            var_0:_ndim = 3 ;
 *            var_0:_dimids = 0, 1, 2 ;
 *            var_0:_datatype = 4 ;
 *            var_0:_varkind = 1 ;
 *            var_0:_chunkdim = 1, 5, 5 ;
 *            var_0:_filter = 2 ;
 *            var_0:_metaoffset = 8LL ;
 *        int var_1 ;
 *            var_1:_ndim = 3 ;
 *            var_1:_dimids = 0, 1, 2 ;
 *            var_1:_datatype = 4 ;
 *            var_1:_varkind = 1 ;
 *            var_1:_chunkdim = 1, 5, 5 ;
 *            var_1:_filter = 2 ;
 *            var_1:_metaoffset = 65544LL ;
 *        byte _datablock_0(_datablock_dim_0) ;
 *            _datablock_0:_varkind = 2 ;
 *        byte _datablock_1(_datablock_dim_1) ;
 *            _datablock_1:_varkind = 2 ;
 *
 *    // global attributes:
 *            :_comressed = 1 ;
 *            :_nwrite = 2 ;
 *            :_recsize = 2LL ;
 *    }
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <mpi.h>
#include <pnetcdf.h>

#define NTIMES 2
#define NY 10
#define NX 10
#define NVARS 2

static int verbose;

#define PNC_ERR(fname) { \
    if (err != NC_NOERR) { \
        printf("Error at %s:%d when calling %s (%s)\n", __FILE__,__LINE__, \
        fname, ncmpi_strerror(err)); \
        nerrs++; \
        goto err_out; \
    } \
}

#define MPI_ERROR(fname) { \
    if (err != MPI_SUCCESS) { \
        int errorStringLen; \
        char errorString[MPI_MAX_ERROR_STRING]; \
        MPI_Error_string(err, errorString, &errorStringLen); \
        printf("Error at %s:%d when calling %s (%s)\n", __FILE__,__LINE__, \
        fname, errorString); \
        nerrs++; \
        goto err_out; \
    } \
}

#define CALC_START_COUNT(len, nprocs, rank, start, count) { \
    count = len / nprocs; \
    start = count * rank; \
    if (rank < len % nprocs) { \
        start += rank; \
        count++; \
    } \
    else { \
        start += len % nprocs; \
    } \
}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] [-k format] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [-k format] file format: 1 for CDF-1, 2 for CDF-2, 5 for CDF-5\n"
    "       [filename] output netCDF file name\n";
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

/*----< compress() >--------------------------------------------------------*/
static int
compress(MPI_Comm comm, char *filename, int cmode)
{
    char name[64];
    int i, j, rank, nprocs, err, nerrs=0, ncid, varid[NVARS];
    int dimid[3], psize[2], rank_y, rank_x;
    MPI_Offset  global_ny, global_nx;
    MPI_Offset start[3], count[3];
    MPI_Info info;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    /* Creates a division of processors in a Cartesian grid */
    psize[0] = psize[1] = 0;
    err = MPI_Dims_create(nprocs, 2, psize);
    MPI_ERROR("MPI_Dims_create");
    if (verbose && rank == 0)
        printf("MPI_Dims_create() 2D: psize=%d %d\n", psize[0],psize[1]);

    /* set rank along X and Y */
    rank_y = rank / psize[1];
    rank_x = rank % psize[1];
    if (verbose && rank == 0)
        printf("Local rank 2D: rank_y=%d rank_x=%d\n", rank_y, rank_x);

    /* set chunking (1st dimension should always be 1 for record variable) */
    int chunk_dim[3];
    chunk_dim[0] = 1;
    chunk_dim[1] = NY / psize[0];
    if (NY % psize[0]) chunk_dim[1]++;
    chunk_dim[2] = NX / psize[1];
    if (NX % psize[1]) chunk_dim[2]++;
    if (verbose && rank == 0)
        printf("chunk_dim: %d %d %d\n", chunk_dim[0],chunk_dim[1],chunk_dim[2]);

    /* set subarray start and count. Each rank writes a subarray of size
     * count[0] x count[1] from offset start[0], start[1], a checkerboard
     * partitioning pattern.
     */
    CALC_START_COUNT(NY, psize[0], rank_y, start[1], count[1])
    CALC_START_COUNT(NY, psize[1], rank_x, start[2], count[2])
    start[0] = 0;
    count[0] = 1;
    if (verbose)
        printf("rank %d: start=%lld %lld %lld count=%lld %lld %lld\n", rank,
               start[0],start[1],start[2], count[0],count[1],count[2]);

    /* allocate write buffer of size count[1] x count[2] */
    int *buf = (int*) malloc(sizeof(int) * count[1] * count[2]);
    for (i=0; i<count[1]*count[2]; i++) buf[i] = rank + i;

    /* create a new file and enable chunking and compression ----------------*/
    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_chunking", "enable");
    MPI_Info_set(info, "nc_chunk_default_filter", "zlib");

    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, info, &ncid);
    PNC_ERR("ncmpi_create")

    /* define dimensions x and y */
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0]);
    PNC_ERR("ncmpi_def_dim")
    err = ncmpi_def_dim(ncid, "Y",    NY,           &dimid[1]);
    PNC_ERR("ncmpi_def_dim")
    err = ncmpi_def_dim(ncid, "X",    NX,           &dimid[2]);
    PNC_ERR("ncmpi_def_dim")

    /* define 3D variables of integer type */
    for (i=0; i<NVARS; i++) {
        snprintf(name, 64, "var_%d", i);
        err = ncmpi_def_var(ncid, name, NC_INT, 3, dimid, &varid[i]);
        PNC_ERR("ncmpi_def_var")

        /* set chunking */
        err = ncmpi_var_set_chunk(ncid, varid[i], chunk_dim);;
        PNC_ERR("ncmpi_var_set_chunk")

        /* use default filter */
        err = ncmpi_var_set_filter(ncid, varid[i], NC_FILTER_DEFLATE);
        PNC_ERR("ncmpi_var_set_filter")
    }

    /* exit define mode */
    err = ncmpi_enddef(ncid);
    PNC_ERR("ncmpi_enddef")

    /* write one time record at a time */
    for (i=0; i<NTIMES; i++) {
        start[0] = i;
        for (j=0; j<NVARS; j++) {
            err = ncmpi_iput_vara_int(ncid, varid[j], start, count, buf, NULL);
            PNC_ERR("ncmpi_iput_vara_int")
        }

        /* wait for nonblocking request to complete */
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
        PNC_ERR("ncmpi_wait_all")
    }

    /* check the current record dimension size */
    MPI_Offset dim_len;
    err = ncmpi_inq_dimlen(ncid, 0, &dim_len);
    PNC_ERR("ncmpi_inq_dimlen")
    if (verbose && rank == 0)
        printf("Time dimension length (expect %lld) and got %lld\n",
                NTIMES * nprocs, dim_len);

    err = ncmpi_close(ncid);
    PNC_ERR("ncmpi_close")

    MPI_Info_free(&info);

    free(buf);

err_out:
    return nerrs;
}

/*----< decompress() >-------------------------------------------------------*/
static int
decompress(MPI_Comm comm, char *filename)
{
    char name[64];
    int i, j, rank, nprocs, err, nerrs=0, ncid, *varid, ulimit_dimid;
    int nvars, dimids[3], filter, chunk_dim[3], psize[2], rank_y, rank_x;
    MPI_Offset nrecs, global_ny, global_nx;
    MPI_Offset start[3], count[3];
    MPI_Info info;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    /* Creates a division of processors in a Cartesian grid */
    psize[0] = psize[1] = 0;
    err = MPI_Dims_create(nprocs, 2, psize);
    MPI_ERROR("MPI_Dims_create");
    if (verbose && rank == 0)
        printf("MPI_Dims_create() 2D: psize=%d %d\n", psize[0],psize[1]);

    /* set rank along X and Y */
    rank_y = rank / psize[1];
    rank_x = rank % psize[1];
    if (verbose && rank == 0)
        printf("Local rank 2D: rank_y=%d rank_x=%d\n", rank_y, rank_x);

    /* open the file for reading with chunking and compression enabled */
    MPI_Info_create(&info);
    MPI_Info_set(info, "nc_chunking", "enable");

    err = ncmpi_open(comm, filename, NC_NOWRITE, info, &ncid);
    PNC_ERR("ncmpi_open")

    MPI_Info_free(&info);

    /* obtain dimension info */
    err = ncmpi_inq_dimid(ncid, "Y", &dimids[1]);
    PNC_ERR("ncmpi_inq_dimlen")

    err = ncmpi_inq_dimid(ncid, "X", &dimids[2]);
    PNC_ERR("ncmpi_inq_dimlen")

    err = ncmpi_inq_dimlen(ncid, dimids[1], &global_ny);
    PNC_ERR("ncmpi_inq_dimlen")

    err = ncmpi_inq_dimlen(ncid, dimids[2], &global_nx);
    PNC_ERR("ncmpi_inq_dimlen")

    /* obtain the number of record variables */
    err = ncmpi_inq_num_rec_vars(ncid, &nvars);
    PNC_ERR("ncmpi_inq_num_rec_vars")
    if (verbose && rank == 0)
        printf("Number of record variables = %d\n", nvars);

    varid = (int*) malloc(sizeof(int) * nvars);

    /* obtain variable ID and dimension info */
    for (i=0; i<nvars; i++) {
        snprintf(name, 64, "var_%d", i);
        err = ncmpi_inq_varid(ncid, name, &varid[i]);
        PNC_ERR("ncmpi_inq_varid")
    }

    /* obtain unlimited dimension ID */
    err = ncmpi_inq_unlimdim(ncid, &ulimit_dimid);
    PNC_ERR("ncmpi_inq_unlimdim")

    /* check the current record dimension size */
    err = ncmpi_inq_dimlen(ncid, ulimit_dimid, &nrecs);
    PNC_ERR("ncmpi_inq_dimlen")
    if (verbose && rank == 0)
        printf("Time dimension length = %lld\n", nrecs);

    /* get chunking */
    for (i=0; i<nvars; i++) {
        err = ncmpi_inq_varname(ncid, varid[i], name);;
        PNC_ERR("ncmpi_inq_varname")
        err = ncmpi_var_get_chunk(ncid, varid[i], chunk_dim);;
        PNC_ERR("ncmpi_var_get_chunk")
        if (verbose && rank == 0)
            printf("var %s chunk_dim[3]=%d %d %d\n", name,
                chunk_dim[0],chunk_dim[1],chunk_dim[2]);

        /* get filter */
        err = ncmpi_var_get_filter(ncid, varid[i], &filter);
        PNC_ERR("ncmpi_var_get_filter")
        if (verbose && rank == 0)
            printf("var %s filter is %s\n", name,
                (filter == NC_FILTER_DEFLATE) ?
                "NC_FILTER_DEFLATE": (filter == NC_FILTER_SZ) ?
                "NC_FILTER_SZ" : "UNKNOWN");
    }

    /* set subarray start and count. Each rank reads a subarray of size
     * count[0] x count[1] from offset start[0], start[1], a checkerboard
     * partitioning pattern.
     */
    CALC_START_COUNT(global_ny, psize[0], rank_y, start[1], count[1])
    CALC_START_COUNT(global_nx, psize[1], rank_x, start[2], count[2])
    start[0] = 0;
    count[0] = 1;
    if (verbose)
        printf("rank %d: start=%lld %lld %lld count=%lld %lld %lld\n", rank,
               start[0],start[1],start[2], count[0],count[1],count[2]);

    /* allocate read buffers, one for each variable */
    size_t buf_size = nrecs * count[1] * count[2];
    int **buf = (int**) malloc(sizeof(int*) * nvars);
    for (i=0; i<nvars; i++)
        buf[i] = (int*) malloc(sizeof(int) * buf_size);

    /* Each rank reads a subarray, one at a time. */
    for (i=0; i<nrecs; i++) {
        start[0] = i;
        for (j=0; j<nvars; j++) {
            int *ptr = buf[j] + i * count[1] * count[2];
            err = ncmpi_iget_vara_int(ncid, varid[j], start, count, ptr, NULL);
            PNC_ERR("ncmpi_iget_vara_int")
        }
    }

    /* wait for nonblocking request to complete */
    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
    PNC_ERR("ncmpi_wait_all")

    err = ncmpi_close(ncid);
    PNC_ERR("ncmpi_close")

    for (i=0; i<nvars; i++) free(buf[i]);
    free(buf);
    free(varid);

err_out:
    return nerrs;
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    char filename[256];
    int i, rank, kind=0, cmode=0, err;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 1;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hqk:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'k': kind = atoi(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (verbose && rank == 0) printf("%s: example of using compression\n",__FILE__);

    switch (kind) {
        case(2): cmode = NC_64BIT_OFFSET; break;
        case(5): cmode = NC_64BIT_DATA;   break;
        default: cmode = 0;
    }

    err = compress(MPI_COMM_WORLD, filename, cmode);
    if (err != 0) goto err_out;

    err = decompress(MPI_COMM_WORLD, filename);
    if (err != 0) goto err_out;

    err = pnetcdf_check_mem_usage(MPI_COMM_WORLD);
    if (err != 0) goto err_out;

err_out:
    MPI_Finalize();
    return (err != 0);
}

