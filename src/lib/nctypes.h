/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifndef NCTYPES_H
#define NCTYPES_H

#ifndef HAVE_PTRDIFF_T
typedef int ptrdiff_t;
#endif

#ifndef HAVE_SSIZE_T
typedef int ssize_t;
#endif

#ifndef HAVE_UCHAR
typedef unsigned char uchar;
#endif

#ifndef HAVE_USHORT
typedef unsigned short int  ushort;
#endif

#ifndef HAVE_UINT
typedef unsigned       int  uint;
#endif

#ifndef HAVE_INT64
typedef          long long  int64;
#endif

#ifndef HAVE_UINT64
typedef unsigned long long  uint64;
#endif

#endif
