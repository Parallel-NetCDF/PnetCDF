/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */
#ifndef _RND_H
#define _RND_H

/* useful for aligning memory */
#define	_RNDUP(x, unit)		((((x) + (unit) - 1) / (unit)) * (unit))
#define	_RNDDOWN(x, unit)	((x) - ((x)%(unit)))

/* #define M_RND_UNIT	(sizeof(double))
 * SIZEOF_DOUBLE is defined in config.h
 */
#define M_RND_UNIT	SIZEOF_DOUBLE
#define	M_RNDUP(x) 	_RNDUP(x, M_RND_UNIT)
#define	M_RNDDOWN(x)	__RNDDOWN(x, M_RND_UNIT)

#define D_RNDUP(x, align) _RNDUP(x, (off_t)(align))

#endif
