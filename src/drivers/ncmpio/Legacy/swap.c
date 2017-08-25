/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <mpi.h>

#include "nc.h"
#include "ncx.h"

#ifndef WORDS_BIGENDIAN
/* LITTLE_ENDIAN: DEC and intel */
/*
 * Routines to convert to BIG ENDIAN.
 * Optimize the swapn?b() and swap?b() routines aggressively.
 */

/* netcdf-3.6.2beta5 added loop unrolling to many of these routines.  Could
 * confer a 22% performance increase on little endian platforms if compiler
 * does not already aggressively unroll loops */

void
ncmpii_swapn2b(void *dst, const void *src, MPI_Offset nn)
{
	char *op = (char*) dst;
	const char *ip = (char*) src;

	while(nn > 3)
	{
		*op++ = *(++ip);
		*op++ = *(ip++ -1);
		*op++ = *(++ip);
		*op++ = *(ip++ -1);
		*op++ = *(++ip);
		*op++ = *(ip++ -1);
		*op++ = *(++ip);
		*op++ = *(ip++ -1);
		nn -= 4;
	}
	while(nn-- != 0)
	{
		*op++ = *(++ip);
		*op++ = *(ip++ -1);
	}
}

# ifndef vax
void
ncmpii_swap4b(void *dst, const void *src)
{
	char *op = (char*) dst;
	const char *ip = (char*) src;
	op[0] = ip[3];
	op[1] = ip[2];
	op[2] = ip[1];
	op[3] = ip[0];
}
# endif /* !vax */

void
ncmpii_swapn4b(void *dst, const void *src, MPI_Offset nn)
{
	char *op = (char*) dst;
	const char *ip = (char*) src;

	while(nn > 3)
	{
		op[0] = ip[3];
		op[1] = ip[2];
		op[2] = ip[1];
		op[3] = ip[0];
		op[4] = ip[7];
		op[5] = ip[6];
		op[6] = ip[5];
		op[7] = ip[4];
		op[8] = ip[11];
		op[9] = ip[10];
		op[10] = ip[9];
		op[11] = ip[8];
		op[12] = ip[15];
		op[13] = ip[14];
		op[14] = ip[13];
		op[15] = ip[12];
		op += 16;
		ip += 16;
		nn -= 4;
	}
	while(nn-- != 0)
	{
		op[0] = ip[3];
		op[1] = ip[2];
		op[2] = ip[1];
		op[3] = ip[0];
		op += 4;
		ip += 4;
	}
}

# ifndef vax
void
ncmpii_swap8b(void *dst, const void *src)
{
	char *op = (char*) dst;
	const char *ip = (char*) src;
	op[0] = ip[7];
	op[1] = ip[6];
	op[2] = ip[5];
	op[3] = ip[4];
	op[4] = ip[3];
	op[5] = ip[2];
	op[6] = ip[1];
	op[7] = ip[0];
}
# endif /* !vax */

# ifndef vax
void
ncmpii_swapn8b(void *dst, const void *src, MPI_Offset nn)
{
	char *op = (char*) dst;
	const char *ip = (char*) src;

	while(nn > 1)
	{
		op[0] = ip[7];
		op[1] = ip[6];
		op[2] = ip[5];
		op[3] = ip[4];
		op[4] = ip[3];
		op[5] = ip[2];
		op[6] = ip[1];
		op[7] = ip[0];
		op[8] = ip[15];
		op[9] = ip[14];
		op[10] = ip[13];
		op[11] = ip[12];
		op[12] = ip[11];
		op[13] = ip[10];
		op[14] = ip[9];
		op[15] = ip[8];
		op += 16;
		ip += 16;
		nn -= 2;
	}
	while(nn-- != 0)
	{
		op[0] = ip[7];
		op[1] = ip[6];
		op[2] = ip[5];
		op[3] = ip[4];
		op[4] = ip[3];
		op[5] = ip[2];
		op[6] = ip[1];
		op[7] = ip[0];
		op += 8;
		ip += 8;
	}
}
# endif /* !vax */

#endif /* LITTLE_ENDIAN */

