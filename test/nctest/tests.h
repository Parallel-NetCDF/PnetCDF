/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header$
 *********************************************************************/

#undef PROTO
#ifndef NO_HAVE_PROTOTYPES 
#   define	PROTO(x)	x
#else
#   define	PROTO(x)	()
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern void	test_nccreate	PROTO((
				       char*
				       ));
extern void	test_ncopen	PROTO((
				       char*
				       ));
extern void	test_ncredef	PROTO((
				       char*
				       ));
extern void	test_ncendef	PROTO((
				       char*
				       ));
extern void	test_ncclose	PROTO((
				       char*
				       ));
extern void	test_ncinquire	PROTO((
				       char*
				       ));
extern void	test_ncsync	PROTO((
				       char*
				       ));
extern void	test_ncabort	PROTO((
				       char*
				       ));
extern void	test_ncdimdef	PROTO((
				       char*
				       ));
extern void	test_ncdimid	PROTO((
				       char*
				       ));
extern void	test_ncdiminq	PROTO((
				       char*
				       ));
extern void	test_ncdimrename	PROTO((
					       char*
					       ));
extern void	test_ncvardef	PROTO((
				       char*
				       ));
extern void	test_ncvarid	PROTO((
				       char*
				       ));
extern void	test_ncvarinq	PROTO((
				       char*
				       ));
extern void	test_ncvarput1	PROTO((
				       char*
				       ));
extern void	test_ncvarget1	PROTO((
				       char*
				       ));
extern void	test_ncvarput	PROTO((
				       char*
				       ));
extern void	test_ncvarget	PROTO((
				       char*
				       ));
extern void	test_ncvarputg	PROTO((
				       char*
				       ));
extern void	test_ncvargetg	PROTO((
				       char*
				       ));
extern void	test_ncrecinq	PROTO((
				       char*
				       ));
extern void	test_ncrecput	PROTO((
				       char*
				       ));
extern void	test_ncrecget	PROTO((
				       char*
				       ));
extern void	test_ncvarrename	PROTO((
					       char*
					       ));
extern void	test_ncattput	PROTO((
				       char*
				       ));
extern void	test_ncattinq	PROTO((
				       char*
				       ));
extern void	test_ncattget	PROTO((
				       char*
				       ));
extern void	test_ncattcopy	PROTO((
				       char*,
				       char*
				       ));
extern void	test_ncattname	PROTO((
				       char*
				       ));
extern void	test_ncattrename	PROTO((
					       char*
					       ));
extern void	test_ncattdel	PROTO((
				       char*
				       ));
extern void	test_nctypelen	PROTO((
					void
				       ));
extern int	test_varputget	PROTO((
				       int
				       ));
extern int	test_varputgetg	PROTO((
				       int
				       ));
extern int	test_slabs	PROTO((
				       int
				       ));

#ifdef __cplusplus
}
#endif
