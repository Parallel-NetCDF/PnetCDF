/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#ifndef HAVE_PTRDIFF_T
typedef long int ptrdiff_t;
#endif

#ifndef HAVE_SCHAR
typedef signed char schar;
#endif

#ifndef HAVE_UCHAR
typedef unsigned char uchar;
#endif

#ifndef HAVE_USHORT
typedef unsigned short int ushort;
#endif

#ifndef HAVE_UINT
typedef unsigned int uint;
#endif

#ifndef HAVE_LONGLONG
typedef long long longlong;
#endif

#ifndef HAVE_ULONGLONG
typedef unsigned long long ulonglong;
#endif

#ifndef HAVE_INT64
typedef long long int64;
#endif

#ifndef HAVE_UINT64
typedef unsigned long long uint64;
#endif

