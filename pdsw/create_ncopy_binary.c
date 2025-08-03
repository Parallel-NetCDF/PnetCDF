/*********************************************************************
 *
 *  Copyright (C) 2023, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <pnetcdf.h>
#include <assert.h>
#include "baseline_ncx_app.h"

#define NCOPY 10

static int verbose;
const char *source_file = NULL;
const char *output_name = NULL;

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}


int build_ncopy_hdr(struct hdr *input_hdr, struct hdr *output_hdr) {
    int nerrs = 0;
    int ndims = input_hdr->dims.ndefined;
    int nvars = input_hdr->vars.ndefined;

    output_hdr->dims.ndefined = NCOPY * ndims;
    output_hdr->vars.ndefined = NCOPY * nvars;
    output_hdr->dims.value = (hdr_dim **) malloc(output_hdr->dims.ndefined * sizeof(hdr_dim *));
    output_hdr->vars.value = (hdr_var **) malloc(output_hdr->vars.ndefined * sizeof(hdr_var *));
    output_hdr->xsz = 4 * sizeof(uint32_t);  // Tags + counts for dims and vars

    // ----------------- Duplicate Dimensions -----------------
    for (int n = 0; n < NCOPY; n++) {
        for (int i = 0; i < ndims; i++) {
            int index = n * ndims + i;
            output_hdr->dims.value[index] = (hdr_dim *) malloc(sizeof(hdr_dim));
            hdr_dim *src = input_hdr->dims.value[i];
            hdr_dim *dst = output_hdr->dims.value[index];

            char name[256];
            snprintf(name, sizeof(name), "%s_copy_%d", src ->name, n);
            dst->name = strdup(name);
            dst->name_len = strlen(dst->name);
            dst->size = src->size;

            output_hdr->xsz += sizeof(uint32_t) + dst->name_len + sizeof(uint32_t);
        }
    }

    // ----------------- Duplicate Variables -----------------
    for (int n = 0; n < NCOPY; n++) {
        for (int i = 0; i < nvars; i++) {
            int index = n * nvars + i;
            output_hdr->vars.value[index] = (hdr_var *) malloc(sizeof(hdr_var));
            hdr_var *src = input_hdr->vars.value[i];
            hdr_var *dst = output_hdr->vars.value[index];

            char name[256];
            snprintf(name, sizeof(name), "copy_%d_%s", n + 1, src->name);
            dst->name = strdup(name);
            dst->name_len = strlen(dst->name);
            dst->xtype = src->xtype;
            dst->ndims = src->ndims;

            dst->dimids = (int *) malloc(dst->ndims * sizeof(int));
            for (int j = 0; j < dst->ndims; j++) {
                dst->dimids[j] = n * ndims + src->dimids[j];
            }

            dst->attrs.ndefined = src->attrs.ndefined;
            dst->attrs.value = (hdr_attr **) malloc(dst->attrs.ndefined * sizeof(hdr_attr *));

            for (int k = 0; k < dst->attrs.ndefined; k++) {
                hdr_attr *src_attr = src->attrs.value[k];
                hdr_attr *dst_attr = (hdr_attr *) malloc(sizeof(hdr_attr));

                dst_attr->name = strdup(src_attr->name);
                dst_attr->name_len = strlen(dst_attr->name);
                dst_attr->xtype = src_attr->xtype;
                dst_attr->nelems = src_attr->nelems;

                int elem_size;
                xlen_nc_type_meta(dst_attr->xtype, &elem_size);
                dst_attr->xvalue = malloc(elem_size * dst_attr->nelems);
                memcpy(dst_attr->xvalue, src_attr->xvalue, elem_size * dst_attr->nelems);
                dst->attrs.value[k] = dst_attr;

                output_hdr->xsz += sizeof(uint32_t) + dst_attr->name_len;
                output_hdr->xsz += 2 * sizeof(uint32_t);
                output_hdr->xsz += elem_size * dst_attr->nelems;
            }

            output_hdr->xsz += sizeof(uint32_t) + dst->name_len;
            output_hdr->xsz += sizeof(uint32_t); // xtype
            output_hdr->xsz += sizeof(uint32_t); // ndims
            output_hdr->xsz += dst->ndims * sizeof(uint32_t);
            output_hdr->xsz += 2 * sizeof(uint32_t); // tag + attr count
        }
    }

    return nerrs;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) printf("Reading and duplicating metadata...\n");
    if (argc < 3) {
            if (rank == 0)
                fprintf(stderr, "Usage: %s <source_file> <output_file>\n", argv[0]);
            MPI_Finalize();
            return 1;
        }
        source_file = argv[1];
        output_name = argv[2];

    struct hdr input_hdr;
    struct hdr output_hdr;

    read_metadata_from_file(source_file, &input_hdr);
    build_ncopy_hdr(&input_hdr, &output_hdr);

    // ------------------ Serialize and Write ------------------
    char* send_buffer = (char*) malloc(output_hdr.xsz);
    int status = serialize_hdr_meta(&output_hdr, send_buffer);

    if (rank == 0) {
        printf("Output buffer size: %lld bytes\n", output_hdr.xsz);
        FILE* file = fopen(output_name, "wb");
        if (file == NULL) {
            perror("Failed to open output file");
            goto fn_exit;
        }
        size_t written = fwrite(send_buffer, 1, output_hdr.xsz, file);
        if (written != output_hdr.xsz) {
            perror("Failed to write complete buffer to file");
        }
        fclose(file);
    }

fn_exit:
    free_hdr_meta(&input_hdr);
    free_hdr_meta(&output_hdr);
    free(send_buffer);
    MPI_Finalize();
    return 0;
}