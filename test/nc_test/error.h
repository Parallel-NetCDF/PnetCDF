/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */


#ifdef __cplusplus
extern "C" {
#endif

/* Print error message to stderr, don't exit */
extern void	error (const char *fmt, ...)
#ifdef _GNUC_
__attribute__ ((format (printf, 1, 2)))
#endif
;


void print(const char *fmt, ...)
#ifdef _GNUC_
__attribute__ ((format (printf, 1, 2)))
#endif
;


extern int ifFail(const int expr, const int line, const char *file, const char *func);

extern void
print_n_size_t(int nelems, const MPI_Offset *array);

#ifdef __cplusplus
}
#endif

#define IF(EXPR) if (ifFail(EXPR, __LINE__, __FILE__, __func__))
#define ELSE_NOK else {nok++;}
