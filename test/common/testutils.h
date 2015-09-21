/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */


#ifndef _UTILS_H
#define _UTILS_H

#include <limits.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

typedef struct {
	char infname[PATH_MAX];
	char outfname[PATH_MAX];
} params;

void parse_read_args(int argc, char **argv, int rank, params *p);
void parse_write_args(int argc, char **argv, int rank, params *p);

extern char* nc_err_code_name(int err);

#endif
