AC_REVISION($Id$)dnl
dnl -*- Mode: shell-script-mode; -*-
dnl Process this file with GNU autoconf(1) to produce a configure script.
dnl

AC_PREREQ(2.69)
AC_INIT([PnetCDF-FLASH-IO], [1.0.0], [parallel-netcdf@mcs.anl.gov])

AC_CONFIG_AUX_DIR([./scripts])

AM_INIT_AUTOMAKE([foreign])
dnl AM_MAINTAINER_MODE
AC_CONFIG_MACRO_DIRS([m4])

dnl AM_EXTRA_RECURSIVE_TARGETS macro was introduced into automake 1.13
m4_ifdef([AM_EXTRA_RECURSIVE_TARGETS], [AM_EXTRA_RECURSIVE_TARGETS([tests])])

dnl line below indicates this is not built as part of PnetCDF testsets
AM_CONDITIONAL(BUILD_TESTSETS, [false])

dnl user defined macro for printing messages for debugging
_DEBUG=no
AC_DEFUN([UD_MSG_DEBUG],
         [if test "x${_DEBUG}" = xyes ; then
             AC_MSG_NOTICE(DEBUG: $1)
          fi
         ]
)

dnl an option to use a customized rm command
AC_ARG_VAR(RM, Command for deleting files or directories. @<:@default: rm@:>@)
if test "x${RM}" != x ; then
   AC_MSG_CHECKING(rm )
   dnl ${RM} set by user may contain flags
   RM_CMD=`echo ${RM} | cut -d' ' -f1`
   if ! test -f ${RM_CMD} ; then
      AC_CHECK_PROG([rm_cmd], [${RM_CMD}], [yes], [no])
      if test "x${rm_cmd}" = xno ; then
         RM=rm
      fi
   fi
   AC_MSG_RESULT(using $RM)
else
   RM=rm
fi

AC_ARG_ENABLE(echo,
   [AS_HELP_STRING([--enable-echo],
                   [Turn on strong echoing. @<:@default: no@:>@])],
   [set -x]
)

AC_ARG_VAR(MPIF90, MPI Fortran 90 compiler)
# Check env variable $MPIF90 set by the users
if test "x${MPIF90}" != x ; then
    UD_MSG_DEBUG(MPIF90 is set to $MPIF90)
    if test ! -x "${MPIF90}" ; then
        UD_MSG_DEBUG(MPIF90 is defined as ${MPIF90} but file does not exist)
        # MPIF90 may not be an absolute path, check MPIF90 from $PATH
        MPIF90_ISVALID=""
        AC_PATH_PROGS(MPIF90_ISVALID, [${MPIF90}])
        if test "x${MPIF90_ISVALID}" = x ; then
            AC_MSG_WARN(File ${MPIF90} does not exist)
            MPIF90=""
        else
            MPIF90=$MPIF90_ISVALID
        fi
    fi
fi

dnl if MPIF90 is still empty, search under user's PATH
CANDIDATE_MPIF90="mpif90 mpixlf90_r mpixlf90 blrts_xlf90 mpxlf90_r mpxlf90 mpxlf95_r mpxlf95 ftn mpiifort mpiifc mpifrtpx"
if test "x${MPIF90}" = x ; then
   UD_MSG_DEBUG(find possible MPIF90 under user's PATH)
   AC_PATH_PROGS([MPIF90], [$CANDIDATE_MPIF90])
fi

if test "x${MPIF90}" = x ; then
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
    No MPI F90 compiler can be found. FLASH I/O kernel requires an MPI F90
    compiler. Please specify the location of one with the MPIF90 environment
    variable
   -----------------------------------------------------------------------])
fi

dnl Check if MPIF90 is a valid MPI compiler
FC=${MPIF90}
AC_PROG_FC

dnl FCFLAGS is set in AC_PROG_FC
ACX_FC_MISMATCH([AS_VAR_APPEND([FCFLAGS], [" $acx_cv_fc_mismatch_flag"])])

dnl find compile flags for handling files with .F .F90 .f90 extensions
dnl steal AC_FC_PP_SRCEXT from autoconf V2.69 to make AC_FC_PP_SRCEXT
AC_FC_PP_SRCEXT([F90])
AC_SUBST(FCFLAGS_F90)

FPPFLAGS=${FPPFLAGS-}
AC_SUBST(FPPFLAGS)

dnl compiler command-line define preprocessor flag, result in FC_DEFINE
AC_FC_PP_DEFINE

AC_LANG_PUSH([Fortran])
AC_FC_SRCEXT([f90])
AC_SUBST(FCFLAGS_f90)

AC_COMPILE_IFELSE([AC_LANG_CALL([],[MPI_Comm_rank])],
                  [valid_mpif90=yes],[valid_mpif90=no]
)
AC_LANG_POP([Fortran])
if test "x${valid_mpif90}" = xno ; then
   AC_MSG_ERROR([
      ------------------------------------------------------------
        Invalid MPI Fortran 90 compiler specified: "${MPIF90}"
        A working MPI compiler is required. Please specify the
        location of one either with the MPIF90 environment
        variable or the --with-mpi configure flag
      ------------------------------------------------------------])
fi
AC_SUBST(MPIF90)
FC=$MPIF90

AC_FC_MODULE_FLAG
dnl AC_FC_MODULE_FLAG defines FC_MODINC
UD_MSG_DEBUG([FC_MODINC=$FC_MODINC])

PNETCDF_INSTALL=""
AC_ARG_WITH(pnetcdf,
    [AS_HELP_STRING([--with-pnetcdf=/path/to/implementation],
                    [installation prefix for PnetCDF implementation])],
    if test "x${withval}" = xyes ; then
        AC_MSG_ERROR(--with-pnetcdf must be given a pathname)
    else
        PNETCDF_INSTALL=${withval}
    fi
)

AC_MSG_CHECKING(PnetCDF Fortran library)
if test "x${PNETCDF_INSTALL}" != x ; then
   PNETCDF_INC=${PNETCDF_INSTALL}/include
   PNETCDF_LIB="-L${PNETCDF_INSTALL}/lib"
fi

have_pnetcdf=no

AC_LANG_PUSH([Fortran])
FCFLAGS_save=${FCFLAGS}
LIBS_save=${LIBS}
if test "x${PNETCDF_INC}" != x ; then
   FCFLAGS="${FCFLAGS} ${FPPFLAGS} ${FC_MODINC}${PNETCDF_INC}"
fi
if test "x${PNETCDF_LIB}" != x ; then
   LDFLAGS="${LDFLAGS} ${PNETCDF_LIB}"
fi
LIBS="-lpnetcdf ${LIBS}"
AC_LINK_IFELSE([AC_LANG_CALL([], [nfmpi_create])],
   [AC_MSG_RESULT(works)
    have_pnetcdf=yes],
   [AC_MSG_RESULT(failed to compile a PnetCDF Fortran program)
    unset FC
   ]
)
FCFLAGS=${FCFLAGS_save}
LIBS=${LIBS_save}
AC_LANG_POP([Fortran])
AC_SUBST(PNETCDF_INC)
AC_SUBST(PNETCDF_LIB)

if test "x${have_pnetcdf}" = xno ; then
    AC_MSG_ERROR([
    ------------------------------------------------------------
     Valid PnetCDF Fortran library and header include file are required
     to build the FLASH I/O kernel. Use option
         --with-pnetcdf=/path/to/implementation
     to specify the location of PnetCDF build.  Stopping ...
     Check 'config.log' for more information.
    ------------------------------------------------------------])
fi

LT_PREREQ([2.2])
LT_INIT

AC_ARG_VAR(TESTMPIRUN, [MPI run command for running ptest, , @<:@default: mpiexec@:>@])
AC_ARG_VAR(TESTOUTDIR, [Output file directory for running check/ptest, @<:@default: ./@:>@])
if test "x${TESTMPIRUN}" = x ; then
   dnl set default to mpiexec
   AC_PATH_PROGS([TESTMPIRUN], [mpiexec mpirun])
fi
if test "x${TESTOUTDIR}" = x ; then
    dnl set default to current directory
    TESTOUTDIR=.
else
   FSTYPE_PREFIX=`echo ${TESTOUTDIR} | cut -d: -f1`
   if test $FSTYPE_PREFIX = $TESTOUTDIR ; then
      # no name prefix end with ':'
      FSTYPE_PREFIX=
   else
      # check if name prefix is one of file system types known to ROMIO
      romio_known_fstypes=(ufs nfs xfs pvfs2 gpfs panfs lustre daos testfs ime quobyte)
      known_fstype=
      for pre in $romio_known_fstypes ; do
          if test "$FSTYPE_PREFIX" = $pre ; then
             known_fstype=$pre
             break
          fi
      done
      if test "x$known_fstype" = x ; then
         FSTYPE_PREFIX=
      else
         FSTYPE_PREFIX+=":"
         # remove name prefix from test output directory name
         TESTOUTDIR=${TESTOUTDIR#"$FSTYPE_PREFIX"}
      fi
   fi
fi
AC_SUBST(TESTMPIRUN)
AC_SUBST(TESTOUTDIR)
AC_SUBST(FSTYPE_PREFIX)

AM_CONDITIONAL(BUILD_BENCHMARKS_IN_PNETCDF, [test yes = no])

AC_CONFIG_FILES(Makefile)
AC_OUTPUT
echo "--------------------------------------------------------------------"
echo \
"
   FLASH I/O kernel with PnetCDF I/O method

   Compilers:
              MPIF90      = ${MPIF90}
              FFLAGS      = ${FFLAGS}
              FCFLAGS     = ${FCFLAGS}
              F90FLAGS    = ${F90FLAGS}
              LDFLAGS     = ${LDFLAGS}
              LIBS        = ${LIBS}

   Now type 'make' to produce executable
---------------------------------------------------------------------"

