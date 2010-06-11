#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <pnetcdf.h>

#define HANDLE_ERROR {                                \
    if (status != NC_NOERR)                           \
        printf("Error at line %d (%s)\n", __LINE__,   \
               ncmpi_strerror(status));               \
}

#define HANDLE_DIFF(str) {                                       \
    int doStop, isDiff = (str[0] == '\0') ? 0 : 1;               \
    MPI_Allreduce(&isDiff, &doStop, 1, MPI_INT, MPI_MAX, comm);  \
    if (doStop) {                                                \
        printf("P%d: diff at %s",rank,str);                      \
        MPI_Finalize();                                          \
        return 1;                                                \
    }                                                            \
}

#define CHECK_GLOBAL_ATT_DIFF(type, func, nctype) {                          \
    int   pos, len = attlen1 * sizeof(type);                                 \
    type *b1 = (type *)malloc(len);                                          \
    type *b2 = (type *)malloc(len);                                          \
    status = func(ncid1, NC_GLOBAL, name1, b1);                              \
    HANDLE_ERROR                                                             \
    status = func(ncid2, NC_GLOBAL, name2, b2);                              \
    HANDLE_ERROR                                                             \
    if ((pos = memcmp(b1, b2, len)) != 0)                                    \
        sprintf(str,"attribute[%d] %s: %s buf1 != buf2 at position %d\n",    \
                i,name1,#nctype,pos);                                        \
    HANDLE_DIFF(str)                                                         \
    free(b1);                                                                \
    free(b2);                                                                \
    break;                                                                   \
}

#define CHECK_VAR_ATT_DIFF(type, func, nctype) {                             \
    int   pos, len = attlen1 * sizeof(type);                                 \
    type *b1 = (type *)malloc(len);                                          \
    type *b2 = (type *)malloc(len);                                          \
    status = func(ncid1, i, name1, b1);                                      \
    HANDLE_ERROR                                                             \
    status = func(ncid2, i, name2, b2);                                      \
    HANDLE_ERROR                                                             \
    if ((pos = memcmp(b1, b2, len)) != 0)                                    \
        sprintf(str,"variable[%d] %s: attribute[%d] %s: %s buf1 != buf2 at position %d\n", \
                i,name,j,name1,#nctype,pos);                                 \
    HANDLE_DIFF(str)                                                         \
    free(b1);                                                                \
    free(b2);                                                                \
    break;                                                                   \
}

#define CHECK_VAR_DIFF(type, func, nctype) {                                 \
    int   pos, len = varsize * sizeof(type);                                 \
    type *b1 = (type *)malloc(len);                                          \
    type *b2 = (type *)malloc(len);                                          \
    status = func(ncid1, i, start, shape, b1);                               \
    HANDLE_ERROR                                                             \
    status = func(ncid2, i, start, shape, b2);                               \
    HANDLE_ERROR                                                             \
    if ((pos = memcmp(b1, b2, len)) != 0)                                    \
        sprintf(str,"variable[%d] %s: %s buf1 != buf2 at position %d\n",     \
                i,name,#nctype,pos);                                        \
    HANDLE_DIFF(str)                                                         \
    free(b1);                                                                \
    free(b2);                                                                \
    break;                                                                   \
}


int main(int argc, char **argv) {
    int i, j, status, isRecvar, rank, nprocs;
    int ncid1, ndims1, nvars1, natts1, unlimdimid1, dimids1[NC_MAX_DIMS];
    int ncid2, ndims2, nvars2, natts2, unlimdimid2, dimids2[NC_MAX_DIMS];
    char str[512], name1[NC_MAX_NAME], name2[NC_MAX_NAME], name[NC_MAX_NAME];
    MPI_Offset shape[NC_MAX_VAR_DIMS], varsize, start[NC_MAX_VAR_DIMS];
    MPI_Offset attlen1, dimlen1, attlen2, dimlen2;
    ncmpi_type type1, type2;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    if (argc != 3) {
        if (rank == 0)
            printf("Usage: %s file1 file2\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    str[0] = '\0';

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
    if (ndims1 != ndims2)
        sprintf(str,"ndims1(%d) != ndims2(%d)\n",ndims1, ndims2);
    HANDLE_DIFF(str)
    if (nvars1 != nvars2)
        sprintf(str,"nvars1(%d) != nvars2(%d)\n",nvars1, nvars2);
    HANDLE_DIFF(str)
    if (natts1 != natts2)
        sprintf(str,"natts1(%d) != natts2(%d)\n",natts1, natts2);
    HANDLE_DIFF(str)

    /* Inquire global attributes, assume CHAR attributes. */
    for (i=0; i<natts1; i++) {
        status = ncmpi_inq_attname(ncid1, NC_GLOBAL, i, name1);
        HANDLE_ERROR
        status = ncmpi_inq_attname(ncid1, NC_GLOBAL, i, name2);
        HANDLE_ERROR
        if (strcmp(name1, name2) != 0)
            sprintf(str,"attribute[%d] name1(%s) != name2(%s)\n",i,name1,name2);
        HANDLE_DIFF(str)

        status = ncmpi_inq_att(ncid1, NC_GLOBAL, name1, &type1, &attlen1);
        HANDLE_ERROR
        status = ncmpi_inq_att(ncid2, NC_GLOBAL, name2, &type2, &attlen2);
        HANDLE_ERROR
        if (type1 != type2)
            sprintf(str,"attribute[%d] %s: type1(%d) != type2(%d)\n",i,name1,type1,type2);
        HANDLE_DIFF(str)
        if (attlen1 != attlen2)
            sprintf(str,"attribute[%d] %s: attlen1(%lld) != attlen2(%lld)\n",i,name1,(long long int)attlen1,(long long int)attlen2);
        HANDLE_DIFF(str)
        switch (type1) {
            case NC_CHAR:   CHECK_GLOBAL_ATT_DIFF(char,   ncmpi_get_att_text,   NC_CHAR)
            case NC_SHORT:  CHECK_GLOBAL_ATT_DIFF(short,  ncmpi_get_att_short,  NC_SHORT)
            case NC_INT:    CHECK_GLOBAL_ATT_DIFF(int,    ncmpi_get_att_int,    NC_INT)
            case NC_FLOAT:  CHECK_GLOBAL_ATT_DIFF(float,  ncmpi_get_att_float,  NC_FLOAT)
            case NC_DOUBLE: CHECK_GLOBAL_ATT_DIFF(double, ncmpi_get_att_double, NC_DOUBLE)
            default: ; /* TODO: handle unexpected types */
        }
    }

    /* Inquire dimension */
    for (i=0; i<ndims1; i++) {
        status = ncmpi_inq_dim(ncid1, i, name1, &dimlen1);
        HANDLE_ERROR
        status = ncmpi_inq_dim(ncid2, i, name2, &dimlen2);
        HANDLE_ERROR
        if (dimlen1 != dimlen2)
		/* cast to quiet warning on 32 bit platforms */
            sprintf(str,"dimension[%d] %s: dimlen1(%lld) != dimlen2(%lld)\n",i,name1,(long long int)dimlen1,(long long int)dimlen2);
        HANDLE_DIFF(str)
    }

    /* Inquire variables */
    for (i=0; i<nvars1; i++) {
        status = ncmpi_inq_var(ncid1, i, name1, &type1, &ndims1, dimids1, &natts1);
        HANDLE_ERROR
        status = ncmpi_inq_var(ncid2, i, name2, &type2, &ndims2, dimids2, &natts2);
        HANDLE_ERROR
        if (strcmp(name1, name2) != 0)
            sprintf(str,"variable[%d]: name1(%s) != name2(%s)\n",i,name1,name2);
        HANDLE_DIFF(str)
        if (type1 != type2)
            sprintf(str,"variable[%d] %s: type1(%d) != type2(%d)\n",i,name1,type1,type2);
        HANDLE_DIFF(str)
        if (ndims1 != ndims2)
            sprintf(str,"variable[%d] %s: ndims1(%d) != ndims2(%d)\n",i,name1,ndims1,ndims2);
        HANDLE_DIFF(str)
        for (j=0; j<ndims1; j++) {
            if (dimids1[j] != dimids2[j])
                sprintf(str,"variable[%d] %s: dimids1[%d]=%d != dimids2[%d]=%d\n",i,name1,j,dimids1[j],j,dimids2[j]);
            HANDLE_DIFF(str)
        }
        if (natts1 != natts2)
            sprintf(str,"variable[%d] %s: natts1(%d) != natts2(%d)\n",i,name1,natts1,natts2);
        HANDLE_DIFF(str)

        strcpy(name,name1);

        /* var attributes, assume CHAR attributes */
        for (j=0; j<natts1; j++) {
            status = ncmpi_inq_attname(ncid1, i, j, name1);
            HANDLE_ERROR
            status = ncmpi_inq_attname(ncid2, i, j, name2);
            HANDLE_ERROR
            if (strcmp(name1, name2) != 0)
                sprintf(str,"variable[%d] %s: attr name[%d] (%s) != (%s)\n",i,name,j,name1,name2);
            HANDLE_DIFF(str)

            status = ncmpi_inq_att(ncid1, i, name1, &type1, &attlen1);
            HANDLE_ERROR
            status = ncmpi_inq_att(ncid2, i, name2, &type2, &attlen2);
            HANDLE_ERROR
            if (type1 != type2)
                sprintf(str,"variable[%d] %s: attr type[%d] (%d) != (%d)\n",i,name,j,type1,type2);
            HANDLE_DIFF(str)
            if (attlen1 != attlen2)
                sprintf(str,"variable[%d] %s: attr attlen[%d] (%lld) != (%lld)\n",i,name,j,(long long int)attlen1,(long long int)attlen2);
            HANDLE_DIFF(str)

            switch (type1) {
                case NC_CHAR:   CHECK_VAR_ATT_DIFF(char,   ncmpi_get_att_text,   NC_CHAR)
                case NC_SHORT:  CHECK_VAR_ATT_DIFF(short,  ncmpi_get_att_short,  NC_SHORT)
                case NC_INT:    CHECK_VAR_ATT_DIFF(int,    ncmpi_get_att_int,    NC_INT)
                case NC_FLOAT:  CHECK_VAR_ATT_DIFF(float,  ncmpi_get_att_float,  NC_FLOAT)
                case NC_DOUBLE: CHECK_VAR_ATT_DIFF(double, ncmpi_get_att_double, NC_DOUBLE)
                default: ; /* TODO: handle unexpected types */
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
        strcpy(name,name1);
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
            case NC_CHAR:   CHECK_VAR_DIFF(char,   ncmpi_get_vara_text_all,   NC_CHAR)
            case NC_SHORT:  CHECK_VAR_DIFF(short,  ncmpi_get_vara_short_all,  NC_SHORT)
            case NC_INT:    CHECK_VAR_DIFF(int,    ncmpi_get_vara_int_all,    NC_INT)
            case NC_FLOAT:  CHECK_VAR_DIFF(float,  ncmpi_get_vara_float_all,  NC_FLOAT)
            case NC_DOUBLE: CHECK_VAR_DIFF(double, ncmpi_get_vara_double_all, NC_DOUBLE)
            default: ; /* TODO: handle unexpected types */
        }
    }

    status = ncmpi_close(ncid1);
    HANDLE_ERROR
    status = ncmpi_close(ncid2);
    HANDLE_ERROR

    if (rank == 0)
        printf("Two files are the same\n");

    MPI_Finalize();
    return 0;
}
