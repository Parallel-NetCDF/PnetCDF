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
			int *nelems,
			int *iscontig_of_ptypes);

int ncmpii_dtype_get_ptype(MPI_Datatype dtype, 
			   MPI_Datatype *ptype,
			   int *iscontig_of_ptypes);

int ncmpii_data_repack(void *inbuf,
                       int incount,
                       MPI_Datatype intype,
                       void *outbuf,
                       int outcount,
                       MPI_Datatype outtype);

#endif
