#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pnetcdf.h>

#define HANDLE_ERROR {                                \
    if (status != NC_NOERR)                           \
        printf("Error at line %s (%s)\n", __LINE__,   \
               ncmpi_strerror(status));               \
}

#define HANDLE_DIFF(str) { \
    printf("Two files are different: %s",str);  \
    MPI_Finalize();  \
    return 1;  \
}

int main(int argc, char **argv) {
    int i, j, status, isRecvar, rank, nprocs, pos;
    int ncid1, ndims1, nvars1, natts1, unlimdimid1, dimids1[NC_MAX_DIMS];
    int ncid2, ndims2, nvars2, natts2, unlimdimid2, dimids2[NC_MAX_DIMS];
    char str[512], name1[NC_MAX_NAME], name2[NC_MAX_NAME];
    MPI_Offset attlen1, dimlen1, shape[NC_MAX_VAR_DIMS], varsize, start[NC_MAX_VAR_DIMS];
    MPI_Offset attlen2, dimlen2;
    ncmpi_type type1, type2;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    if (argc != 3) {
        if (rank == 0)
            printf("Usage: %s file1 file2\n",argv[0]);
        MPI_Finalize();
        return 1;
    }

    status = ncmpi_open(comm, argv[1], NC_NOWRITE, MPI_INFO_NULL, &ncid1);
    HANDLE_ERROR
    status = ncmpi_open(comm, argv[2], NC_NOWRITE, MPI_INFO_NULL, &ncid2);
    HANDLE_ERROR

    /**
     * Inquire the dataset definitions of input dataset AND
     * Add dataset definitions for output dataset.
     */
    status = ncmpi_inq(ncid1, &ndims1, &nvars1, &natts1, &unlimdimid1);
    HANDLE_ERROR
    status = ncmpi_inq(ncid2, &ndims2, &nvars2, &natts2, &unlimdimid2);
    HANDLE_ERROR
    if (ndims1 != ndims2) {
        sprintf(str,"ndims1(%d) != ndims2(%d)\n",ndims1, ndims2);
        HANDLE_DIFF(str)
    }
    if (nvars1 != nvars2) {
        sprintf(str,"nvars1(%d) != nvars2(%d)\n",nvars1, nvars2);
        HANDLE_DIFF(str)
    }
    if (natts1 != natts2) {
        sprintf(str,"natts1(%d) != natts2(%d)\n",natts1, natts2);
        HANDLE_DIFF(str)
    }

    /* Inquire global attributes, assume CHAR attributes. */
    for (i=0; i<natts1; i++) {
        status = ncmpi_inq_attname(ncid1, NC_GLOBAL, i, name1);
        HANDLE_ERROR
        status = ncmpi_inq_attname(ncid1, NC_GLOBAL, i, name2);
        HANDLE_ERROR
        if (strcmp(name1, name2) != 0) {
            sprintf(str,"attribute[%d] name1(%s) != name2(%s)\n",i,name1,name2);
            HANDLE_DIFF(str)
        }

        status = ncmpi_inq_att(ncid1, NC_GLOBAL, name1, &type1, &attlen1);
        HANDLE_ERROR
        status = ncmpi_inq_att(ncid2, NC_GLOBAL, name2, &type2, &attlen2);
        HANDLE_ERROR
        if (type1 != type2) {
            sprintf(str,"attribute[%d] type1(%d) != type2(%d)\n",i,type1,type2);
            HANDLE_DIFF(str)
        }
        if (attlen1 != attlen2) {
            sprintf(str,"attribute[%d] attlen1(%d) != attlen2(%d)\n",i,attlen1,attlen2);
            HANDLE_DIFF(str)
        }

        switch (type1) {
            case NC_CHAR: {
                int   len = attlen1 * sizeof(char);
                char *b1 = (char *)malloc(len);
                char *b2 = (char *)malloc(len);
                status = ncmpi_get_att_text(ncid1, NC_GLOBAL, name1, b1);
                HANDLE_ERROR
                status = ncmpi_get_att_text(ncid2, NC_GLOBAL, name2, b2);
                HANDLE_ERROR
                if ((pos = memcmp(b1, b2, len)) != 0) {
                    sprintf(str,"attribute[%d] NC_CHAR buf1 != buf2 at position %d\n",i,pos);
                    HANDLE_DIFF(str)
                }
                free(b1);
                free(b2);
                break;
            }
            case NC_SHORT: {
                int   len = attlen1 * sizeof(char);
                short *b1 = (short *)malloc(len);
                short *b2 = (short *)malloc(len);
                status = ncmpi_get_att_short(ncid1, NC_GLOBAL, name1, b1);
                HANDLE_ERROR
                status = ncmpi_get_att_short(ncid2, NC_GLOBAL, name2, b2);
                HANDLE_ERROR
                if ((pos = memcmp(b1, b2, len)) != 0) {
                    sprintf(str,"attribute[%d] NC_SHORT buf1 != buf2 at position %d\n",i,pos);
                    HANDLE_DIFF(str)
                }
                free(b1);
                free(b2);
                break;
            }
            case NC_INT: {
                int   len = attlen1 * sizeof(int);
                int *b1 = (int *)malloc(len);
                int *b2 = (int *)malloc(len);
                status = ncmpi_get_att_int(ncid1, NC_GLOBAL, name1, b1);
                HANDLE_ERROR
                status = ncmpi_get_att_int(ncid2, NC_GLOBAL, name2, b2);
                HANDLE_ERROR
                if ((pos = memcmp(b1, b2, len)) != 0) {
                    sprintf(str,"attribute[%d] NC_INT buf1 != buf2 at position %d\n",i,pos);
                    HANDLE_DIFF(str)
                }
                free(b1);
                free(b2);
                break;
            }
            case NC_FLOAT: {
                int   len = attlen1 * sizeof(float);
                float *b1 = (float *)malloc(len);
                float *b2 = (float *)malloc(len);
                status = ncmpi_get_att_float(ncid1, NC_GLOBAL, name1, b1);
                HANDLE_ERROR
                status = ncmpi_get_att_float(ncid2, NC_GLOBAL, name2, b2);
                HANDLE_ERROR
                if ((pos = memcmp(b1, b2, len)) != 0) {
                    sprintf(str,"attribute[%d] NC_INT buf1 != buf2 at position %d\n",i,pos);
                    HANDLE_DIFF(str)
                }
                free(b1);
                free(b2);
                break;
            }
            case NC_DOUBLE: {
                int   len = attlen1 * sizeof(double);
                double *b1 = (double *)malloc(len);
                double *b2 = (double *)malloc(len);
                status = ncmpi_get_att_double(ncid1, NC_GLOBAL, name1, b1);
                HANDLE_ERROR
                status = ncmpi_get_att_double(ncid2, NC_GLOBAL, name2, b2);
                HANDLE_ERROR
                if ((pos = memcmp(b1, b2, len)) != 0) {
                    sprintf(str,"attribute[%d] NC_INT buf1 != buf2 at position %d\n",i,pos);
                    HANDLE_DIFF(str)
                }
                free(b1);
                free(b2);
                break;
            }
            default: ; /* handle unexpected types */
        }
    }

    /* Inquire dimension */
    for (i=0; i<ndims1; i++) {
        status = ncmpi_inq_dim(ncid1, i, name1, &dimlen1);
        HANDLE_ERROR
        status = ncmpi_inq_dim(ncid2, i, name2, &dimlen2);
        HANDLE_ERROR
        if (dimlen1 != dimlen2) {
            sprintf(str,"dimension[%d] dimlen1(%d) != dimlen2(%d)\n",i,dimlen1,dimlen2);
            HANDLE_DIFF(str)
        }
    }

    /* Inquire variables */
    for (i=0; i<nvars1; i++) {
        status = ncmpi_inq_var(ncid1, i, name1, &type1, &ndims1, dimids1, &natts1);
        HANDLE_ERROR
        status = ncmpi_inq_var(ncid2, i, name2, &type2, &ndims2, dimids2, &natts2);
        HANDLE_ERROR
        if (strcmp(name1, name2) != 0) {
            sprintf(str,"variable[%d]: name1(%s) != name2(%s)\n",i,name1,name2);
            HANDLE_DIFF(str)
        }
        if (type1 != type2) {
            sprintf(str,"variable[%d]: type1(%s) != type2(%s)\n",i,type1,type2);
            HANDLE_DIFF(str)
        }
        if (ndims1 != ndims2) {
            sprintf(str,"variable[%d]: ndims1(%s) != ndims2(%s)\n",i,ndims1,ndims2);
            HANDLE_DIFF(str)
        }
        for (j=0; j<ndims1; j++) {
            if (dimids1[j] != dimids2[j]) {
                sprintf(str,"variable[%d]: dimids1[%d]=%d != dimids2[%d]=%d\n",i,j,dimids1[j],j,dimids2[j]);
                HANDLE_DIFF(str)
            }
        }
        if (natts1 != natts2) {
            sprintf(str,"variable[%d]: natts1(%s) != natts2(%s)\n",i,natts1,natts2);
            HANDLE_DIFF(str)
        }

        /* var attributes, assume CHAR attributes */
        for (j=0; j<natts1; j++) {
            status = ncmpi_inq_attname(ncid1, i, j, name1);
            HANDLE_ERROR
            status = ncmpi_inq_attname(ncid2, i, j, name2);
            HANDLE_ERROR
            if (strcmp(name1, name2) != 0) {
                sprintf(str,"variable[%d]: attr name[%d] (%s) != (%s)\n",i,j,name1,name2);
                HANDLE_DIFF(str)
            }

            status = ncmpi_inq_att(ncid1, i, name1, &type1, &attlen1);
            HANDLE_ERROR
            status = ncmpi_inq_att(ncid2, i, name2, &type2, &attlen2);
            HANDLE_ERROR
            if (type1 != type2) {
                sprintf(str,"variable[%d]: attr type[%d] (%s) != (%s)\n",i,j,type1,type2);
                HANDLE_DIFF(str)
            }
            if (attlen1 != attlen2) {
                sprintf(str,"variable[%d]: attr attlen[%d] (%d) != (%d)\n",i,j,attlen1,attlen2);
                HANDLE_DIFF(str)
            }

            switch (type1) {
                case NC_CHAR: {
                    int   len = attlen1 * sizeof(char);
                    char *b1 = (char*)malloc(len);
                    char *b2 = (char*)malloc(len);
                    status = ncmpi_get_att_text(ncid1, i, name1, b1);
                    HANDLE_ERROR
                    status = ncmpi_get_att_text(ncid2, i, name2, b2);
                    HANDLE_ERROR
                    if ((pos = memcmp(b1, b2, len)) != 0) {
                        sprintf(str,"variable[%d] attribute[%d] NC_CHAR buf1 != buf2 at position %d\n",i,j,pos);
                        HANDLE_DIFF(str)
                    }
                    free(b1);
                    free(b2);
                    break;
                }
                case NC_SHORT: {
                    int   len = attlen1 * sizeof(short);
                    short *b1 = (short*)malloc(len);
                    short *b2 = (short*)malloc(len);
                    status = ncmpi_get_att_short(ncid1, i, name1, b1);
                    HANDLE_ERROR
                    status = ncmpi_get_att_short(ncid2, i, name2, b2);
                    HANDLE_ERROR
                    if ((pos = memcmp(b1, b2, len)) != 0) {
                        sprintf(str,"variable[%d] attribute[%d] NC_SHORT buf1 != buf2 at position %d\n",i,j,pos);
                        HANDLE_DIFF(str)
                    }
                    free(b1);
                    free(b2);
                    break;
                }
                case NC_INT: {
                    int   len = attlen1 * sizeof(int);
                    int *b1 = (int*)malloc(len);
                    int *b2 = (int*)malloc(len);
                    status = ncmpi_get_att_int(ncid1, i, name1, b1);
                    HANDLE_ERROR
                    status = ncmpi_get_att_int(ncid2, i, name2, b2);
                    HANDLE_ERROR
                    if ((pos = memcmp(b1, b2, len)) != 0) {
                        sprintf(str,"variable[%d] attribute[%d] NC_INT buf1 != buf2 at position %d\n",i,j,pos);
                        HANDLE_DIFF(str)
                    }
                    free(b1);
                    free(b2);
                    break;
                }
                case NC_FLOAT: {
                    int   len = attlen1 * sizeof(float);
                    float *b1 = (float*)malloc(len);
                    float *b2 = (float*)malloc(len);
                    status = ncmpi_get_att_float(ncid1, i, name1, b1);
                    HANDLE_ERROR
                    status = ncmpi_get_att_float(ncid2, i, name2, b2);
                    HANDLE_ERROR
                    if ((pos = memcmp(b1, b2, len)) != 0) {
                        sprintf(str,"variable[%d] attribute[%d] NC_FLOAT buf1 != buf2 at position %d\n",i,j,pos);
                        HANDLE_DIFF(str)
                    }
                    free(b1);
                    free(b2);
                    break;
                }
                case NC_DOUBLE: {
                    int   len = attlen1 * sizeof(double);
                    double *b1 = (double*)malloc(len);
                    double *b2 = (double*)malloc(len);
                    status = ncmpi_get_att_double(ncid1, i, name1, b1);
                    HANDLE_ERROR
                    status = ncmpi_get_att_double(ncid2, i, name2, b2);
                    HANDLE_ERROR
                    if ((pos = memcmp(b1, b2, len)) != 0) {
                        sprintf(str,"variable[%d] attribute[%d] NC_DOUBLE buf1 != buf2 at position %d\n",i,j,pos);
                        HANDLE_DIFF(str)
                    }
                    free(b1);
                    free(b2);
                    break;
                }
                default:
                  ;
                /* handle unexpected types */
            }
        }
    }

    /**
     * Read data of variables from input dataset 
     * (ONLY DEAL WITH: NC_INT, NC_FLOAT, NC_DOUBLE for now)
     * Write the data out to the corresponding variables in the output dataset
     *
     *  Data Partition (Assume 4 processors):
     *   square: 2-D, (Block, *), 25*100 from 100*100
     *   cube:   3-D, (Block, *, *), 25*100*100 from 100*100*100
     *   xytime: 3-D, (Block, *, *), 25*100*100 from 100*100*100
     *   time:   1-D, Block-wise, 25 from 100
     *
     *  Data Mode API: collective
     */

    for (i=0; i<NC_MAX_VAR_DIMS; i++)
        start[i] = 0;

    for (i=0; i<nvars1; i++) {
        isRecvar = 0;
        varsize = 1;
        ncmpi_inq_var(ncid1, i, name1, &type1, &ndims1, dimids1, &natts1);
        for (j=0; j<ndims1; j++) {
            status = ncmpi_inq_dim(ncid1, dimids1[j], name2, shape + j);
            HANDLE_ERROR
            /* name2 will be discarded */
            if (j == 0) {
                shape[j] /= nprocs;
                start[j] = shape[j] * rank;
            }
            varsize *= shape[j];
            if (dimids1[j] == unlimdimid1)
                isRecvar = 1;
        }
        switch (type1) {
            case NC_CHAR: {
                int   len = varsize * sizeof(char);
                char *b1 = (void *)malloc(len);
                char *b2 = (void *)malloc(len);
                status = ncmpi_get_vara_text_all(ncid1, i, start, shape, b1);
                HANDLE_ERROR
                status = ncmpi_get_vara_text_all(ncid2, i, start, shape, b2);
                HANDLE_ERROR
                if ((pos = memcmp(b1, b2, len) != 0)) {
                    sprintf(str,"variable[%d] %s NC_CHAR buf1 != buf2 at position %d\n",i,name1,pos);
                    HANDLE_DIFF(str)
                }
                free(b1);
                free(b2);
                break; 
            }
            case NC_SHORT: {
                int   len = varsize * sizeof(short);
                short *b1 = (void *)malloc(len);
                short *b2 = (void *)malloc(len);
                status = ncmpi_get_vara_short_all(ncid1, i, start, shape, b1);
                HANDLE_ERROR
                status = ncmpi_get_vara_short_all(ncid2, i, start, shape, b2);
                HANDLE_ERROR
                if ((pos = memcmp(b1, b2, len)) != 0) {
                    sprintf(str,"variable[%d] %s NC_SHORT buf1 != buf2 at position %d\n",i,name1,pos);
                    HANDLE_DIFF(str)
                }
                free(b1);
                free(b2);
                break; 
            }
            case NC_INT: {
                int   len = varsize * sizeof(int);
                int *b1 = (void *)malloc(len);
                int *b2 = (void *)malloc(len);
                status = ncmpi_get_vara_int_all(ncid1, i, start, shape, b1);
                HANDLE_ERROR
                status = ncmpi_get_vara_int_all(ncid2, i, start, shape, b2);
                HANDLE_ERROR
                if ((pos = memcmp(b1, b2, len)) != 0) {
                    sprintf(str,"variable[%d] %s NC_INT buf1 != buf2 at position %d\n",i,name1,pos);
                    HANDLE_DIFF(str)
                }
                free(b1);
                free(b2);
                break; 
            }
            case NC_FLOAT: {
                int   len = varsize * sizeof(float);
                float *b1 = (void *)malloc(len);
                float *b2 = (void *)malloc(len);
                status = ncmpi_get_vara_float_all(ncid1, i, start, shape, b1);
                HANDLE_ERROR
                status = ncmpi_get_vara_float_all(ncid2, i, start, shape, b2);
                HANDLE_ERROR
                if ((pos = memcmp(b1, b2, len)) != 0) {
                    sprintf(str,"variable[%d] %s NC_FLOAT buf1 != buf2 at position %d\n",i,name1,pos);
                    HANDLE_DIFF(str)
                }
                free(b1);
                free(b2);
                break; 
            }
            case NC_DOUBLE: {
                int   len = varsize * sizeof(double);
                double *b1 = (void *)malloc(len);
                double *b2 = (void *)malloc(len);
                status = ncmpi_get_vara_double_all(ncid1, i, start, shape, b1);
                HANDLE_ERROR
                status = ncmpi_get_vara_double_all(ncid2, i, start, shape, b2);
                HANDLE_ERROR
                if ((pos = memcmp(b1, b2, len)) != 0) {
                    sprintf(str,"variable[%d] %s NC_DOUBLE buf1 != buf2 at position %d\n",i,name1,pos);
                    HANDLE_DIFF(str)
                }
                free(b1);
                free(b2);
                break; 
            }
            default:
                ;
                /* handle unexpected types */
        }
    }

    /**
     * Close the datasets
     */
    status = ncmpi_close(ncid1);
    HANDLE_ERROR
    status = ncmpi_close(ncid2);
    HANDLE_ERROR

    if (rank == 0)
        printf("Two files are the same\n");

    MPI_Finalize();
    return 0;
}
