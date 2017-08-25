/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _NCMPIDTYPE_H
#define _NCMPIDTYPE_H

int ncmpio_dtype_decode(MPI_Datatype dtype,
			MPI_Datatype *ptype,
			int *el_size,
			MPI_Offset *nelems,
			int *isderived,
			int *iscontig_of_ptypes);

int ncmpio_data_repack(void *inbuf,
                       MPI_Offset incount,
                       MPI_Datatype intype,
                       void *outbuf,
                       MPI_Offset outcount,
                       MPI_Datatype outtype);

#endif
