/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef _FBITS_H
#define _FBITS_H

/*
 * Macros for dealing with flag bits.
 */
#define fSet(t, f)       ((t) |= (f))
#define fClr(t, f)       ((t) &= ~(f))
#define fIsSet(t, f)     ((t) & (f))
#define fMask(t, f)     ((t) & ~(f))

/*
 * Propositions
 */
/* a implies b */
#define pIf(a,b) (!(a) || (b))
/* a if and only if b, use == when it makes sense */
#define pIff(a,b) (((a) && (b)) || (!(a) && !(b)))

#endif /* _FBITS_H */
