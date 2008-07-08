/*****************************************************************************
 *
 * This file is created by Northwestern University and Argonne National
 * Laboratory
 *
 ****************************************************************************/

#ifndef NCMPI_DTYPE_H
#define NCMPI_DTYPE_H

int ncmpii_dtype_decode(MPI_Datatype dtype,
			MPI_Datatype *ptype,
			int *el_size,
			int64_t *nelems,
			int *isderived,
			int *iscontig_of_ptypes);

int ncmpii_data_repack(void *inbuf,
                       int64_t incount,
                       MPI_Datatype intype,
                       void *outbuf,
                       int64_t outcount,
                       MPI_Datatype outtype);

#endif
