dnl
dnl Copyright (C) 2017, Northwestern University and Argonne National Laboratory
dnl See COPYRIGHT notice in top-level directory.
dnl

dnl -*- Mode: shell-script-mode; -*-
dnl Process this file with GNU autoconf(1) to produce a configure script.
dnl

AC_REVISION([m4_esyscmd_s([git describe --always])])
dnl When using svn, we use SVN keyword Revision to copy the unique revision
dnl number as a stamp into configure script.
dnl AC_REVISION([$Revision$])dnl

dnl autoconf v2.69 was released in 2012-04-24
AC_PREREQ([2.69])
AC_INIT([PnetCDF],[1.12.3],[parallel-netcdf@mcs.anl.gov],[pnetcdf],[https://parallel-netcdf.github.io])

dnl config.h.in will be created by autoreconf (autoheader)
dnl call it right after AC_INIT, as suggested by autoconf
AC_CONFIG_HEADERS([src/include/config.h])
AH_TOP([/*
 * Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 * See COPYRIGHT notice in top-level directory.
 */
])
AC_CONFIG_SRCDIR([src/include/pnetcdf.h.in])

AC_CONFIG_AUX_DIR([./scripts])

dnl call env or printenv to print all user environment variables in config.log
echo '-------------------------------' >& AS_MESSAGE_LOG_FD
echo '  User environment variables:'   >& AS_MESSAGE_LOG_FD
echo '-------------------------------' >& AS_MESSAGE_LOG_FD
env >& AS_MESSAGE_LOG_FD
echo '-------------------------------' >& AS_MESSAGE_LOG_FD
echo '' >& AS_MESSAGE_LOG_FD

dnl Note getting command line should be done before calling AM_INIT_AUTOMAKE
dnl as AM_INIT_AUTOMAKE modifies command line $*
CONFIGURE_ARGS_CLEAN=`echo $* | tr '"' ' '`

dnl AM_INIT_AUTOMAKE([subdir-objects])
AM_INIT_AUTOMAKE([1.13])
dnl enable silent rules by default
AM_SILENT_RULES([yes])

dnl if maintainer mode is disabled, make will *never* attempt to rebuild
dnl configure, Makefile.ins, etc.
AM_MAINTAINER_MODE([enable])

dnl search all local macro files in folder m4
dnl m4_ifdef([AC_CONFIG_MACRO_DIRS], [AC_CONFIG_MACRO_DIRS([m4])], [AC_CONFIG_MACRO_DIR([m4])])
AC_CONFIG_MACRO_DIRS([m4])

dnl AM_EXTRA_RECURSIVE_TARGETS macro was introduced into automake 1.13
dnl m4_ifdef([AM_EXTRA_RECURSIVE_TARGETS], [AM_EXTRA_RECURSIVE_TARGETS([tests])])
AM_EXTRA_RECURSIVE_TARGETS([tests])

dnl parse the version numbers to 4 env variables
PNETCDF_VERSION_MAJOR=`echo ${PACKAGE_VERSION} | cut -d. -f1`
PNETCDF_VERSION_MINOR=`echo ${PACKAGE_VERSION} | cut -d. -f2`
PNETCDF_VERSION_SUB=`echo ${PACKAGE_VERSION} | cut -d. -f3`
PNETCDF_VERSION_PRE=`echo ${PACKAGE_VERSION} | cut -d. -f4`

dnl Note major, minor, and sub are required, but pre is not.
PNETCDF_VERSION=${PACKAGE_VERSION}

dnl Do not change the following line, It is set by SVN automatically.
dnl It defines PNETCDF_RELEASE_DATE, a string that will be used in
dnl ncmpi_inq_libvers() to generate release date
dnl SVN_DATE="$LastChangedDate$"
dnl PNETCDF_RELEASE_DATE2=`echo $SVN_DATE | cut -d' ' -f2`
dnl PNETCDF_RELEASE_DATE=`echo $SVN_DATE | cut -d' ' -f6,7,8 | cut -d')' -f1`

dnl user defined macro for printing messages for debugging
dnl add "pnc_ac_debug=yes" at command line to enable
AC_DEFUN([UD_MSG_DEBUG],
         [if test "x${pnc_ac_debug}" = xyes ; then
             AC_MSG_NOTICE([DEBUG: $1])
          fi
         ]
)

UD_MSG_DEBUG([PNETCDF_VERSION_MAJOR=$PNETCDF_VERSION_MAJOR])
UD_MSG_DEBUG([PNETCDF_VERSION_MINOR=$PNETCDF_VERSION_MINOR])
UD_MSG_DEBUG([PNETCDF_VERSION_SUB=$PNETCDF_VERSION_SUB])
UD_MSG_DEBUG([PNETCDF_VERSION_PRE=$PNETCDF_VERSION_PRE])
UD_MSG_DEBUG([PNETCDF_VERSION=$PNETCDF_VERSION])
dnl UD_MSG_DEBUG([PNETCDF_RELEASE_DATE=$PNETCDF_RELEASE_DATE])

dnl AC_DEFINE_UNQUOTED(PNETCDF_VERSION_MAJOR, $PNETCDF_VERSION_MAJOR, major version number)
dnl AC_DEFINE_UNQUOTED(PNETCDF_VERSION_MINOR, $PNETCDF_VERSION_MINOR, minor version number)
dnl AC_DEFINE_UNQUOTED(PNETCDF_VERSION_SUB, $PNETCDF_VERSION_SUB, sub version number)
dnl AC_DEFINE_UNQUOTED(PNETCDF_VERSION_PRE, $PNETCDF_VERSION_PRE, pre-release string)
dnl AC_DEFINE_UNQUOTED(PNETCDF_VERSION, ["$PNETCDF_VERSION"], full PnetCDF version string)
dnl AC_DEFINE_UNQUOTED(PNETCDF_RELEASE_DATE, ["$PNETCDF_RELEASE_DATE"], PnetCDF release date string)
dnl AC_DEFINE_UNQUOTED(CONFIGURE_ARGS_CLEAN, ["$CONFIGURE_ARGS_CLEAN"], configure command-line arguments used)

AC_SUBST(PNETCDF_VERSION_MAJOR)
AC_SUBST(PNETCDF_VERSION_MINOR)
AC_SUBST(PNETCDF_VERSION_SUB)
AC_SUBST(PNETCDF_VERSION_PRE)
AC_SUBST(PNETCDF_VERSION)

dnl Note that command 'date' is not portable across Unix platforms.
dnl But release date only matters to PnetCDF developers who make the releases.
PNETCDF_RELEASE_DATE="`date '+%B %-d, %Y'`"
AC_SUBST(PNETCDF_RELEASE_DATE)
PNETCDF_RELEASE_DATE_FULL="`date '+%Y-%m-%d'`"
AC_SUBST(PNETCDF_RELEASE_DATE_FULL)

AC_SUBST(PACKAGE_VERSION)
AC_SUBST(CONFIGURE_ARGS_CLEAN)

dnl autoheader only adds these templates to the first invocation of AC_CONFIG_HEADERS
AH_TEMPLATE([ENABLE_FORTRAN],           [Define if to enable Fortran feature])
AH_TEMPLATE([ENABLE_CXX],               [Define if to enable C++ feature])
AH_TEMPLATE([NCBYTE_T],                 [Type of NC_BYTE])
AH_TEMPLATE([NCSHORT_T],                [Type of NC_SHORT])
AH_TEMPLATE([NF_DOUBLEPRECISION_IS_C_], [C type for Fortran double])
AH_TEMPLATE([NF_INT1_IS_C_],            [C type for Fortran INT1])
AH_TEMPLATE([NF_INT1_T],                [Type for Fortran INT1])
AH_TEMPLATE([NF_INT2_IS_C_],            [C type for Fortran INT2])
AH_TEMPLATE([NF_INT2_T],                [Type for Fortran INT2])
AH_TEMPLATE([NF_INT_IS_C_],             [C type for Fortran INT])
AH_TEMPLATE([NF_INT8_IS_C_],            [C type for Fortran INT8])
AH_TEMPLATE([NF_INT8_T],                [Type for Fortran INT8])
AH_TEMPLATE([NF_REAL_IS_C_],            [C type for Fortran REAL])
AH_TEMPLATE([NO_IEEE_FLOAT],            [Does system have IEEE FLOAT])
AH_TEMPLATE([DISABLE_FILE_SYNC],        [Define if to disable MPI_File_sync])
dnl AH_TEMPLATE([ENABLE_IN_PLACE_SWAP], [Define if to enable in-place byte swap])
dnl AH_TEMPLATE([DISABLE_IN_PLACE_SWAP],[Define if to disable in-place byte swap])
AH_TEMPLATE([ENABLE_SUBFILING],         [Define if to enable subfiling feature])
AH_TEMPLATE([ENABLE_NETCDF4],           [Define if to enable NetCDF-4 support])
AH_TEMPLATE([ENABLE_ADIOS],             [Define if to enable ADIOS BP read feature])
AH_TEMPLATE([HDF5_VER_GE_1_10_4],       [Define if HDF5 version is at least 1.10.4])
AH_TEMPLATE([NETCDF_GE_4_5_0],          [Define if NetCDF version is at least 4.5.0])
AH_TEMPLATE([PNC_MALLOC_TRACE],         [Define if to enable malloc tracing])
AH_TEMPLATE([RELAX_COORD_BOUND],        [Define if relaxed coordinate check is enabled])
AH_TEMPLATE([ENABLE_LARGE_SINGLE_REQ],  [Define if to enable large single MPI-IO request])
AH_TEMPLATE([ENABLE_NULL_BYTE_HEADER_PADDING], [Define if to enable strict null-byte padding in file header])
AH_TEMPLATE([ENABLE_BURST_BUFFER],      [Define if to enable burst buffer feature])
AH_TEMPLATE([PNETCDF_PROFILING],        [Define if to enable PnetCDF internal performance profiling])
AH_TEMPLATE([ENABLE_THREAD_SAFE],       [Define if to enable thread-safe capability])
AH_TEMPLATE([ENABLE_REQ_AGGREGATION],   [Define if able to support request aggregation in nonblocking routines])
dnl AH_TEMPLATE([HAVE_MPI_COUNT],       [Define if type MPI_Count is defined])

AH_TOP([#ifndef _CONFIG_H
#define _CONFIG_H])
AH_BOTTOM([#include <nctypes.h>
#endif])

dnl an option to use a customized rm command
RM=rm
AC_ARG_VAR(RM, Command for deleting files or directories. @<:@default: rm@:>@)
if test "x${RM}" != x ; then
   dnl ${RM} set by user may contain flags
   RM_CMD=`echo ${RM} | cut -d' ' -f1`
   # use "if ! test" is not portable, according to autoconf user guide
   if test -f ${RM_CMD} ; then :; else
      AC_CHECK_PROG([rm_cmd], [${RM_CMD}], [yes], [no])
      if test "x${rm_cmd}" = xno ; then
         RM=rm
      fi
   fi
fi
UD_MSG_DEBUG([RM=$RM])

dnl AC_PROG_SED and AC_PROG_GREP are only available on autoconf 2.60 and later
AC_PROG_AWK
AC_PROG_GREP
AC_PROG_EGREP

AC_PROG_SED
dnl check sed command option -i and set SED_I (this requires RM defined)
UD_PROG_SED_I

AC_ARG_ENABLE(echo,
   [AS_HELP_STRING([--enable-echo],
                   [Turn on strong echoing. @<:@default: disabled@:>@])],
   [set -x]
)

MPI_INSTALL=
AC_ARG_WITH(mpi,
   [AS_HELP_STRING([--with-mpi=DIR],
                   [Provide the MPI installation path in DIR.])],
   [ dnl this clause is run when --with-mpi or --without-mpi is used
   if test "x${withval}" = xno ; then
      AC_MSG_ERROR([
      -----------------------------------------------------------------------
        PnetCDF is built on top of MPI. Configure option --without-mpi or
        --with-mpi=no should not be used. Abort.
      -----------------------------------------------------------------------])
   elif test "x${withval}" = x ; then
      # when "--with-mpi=" is used
      AC_MSG_ERROR(--with-mpi is set but the value is NULL)
   elif test "x${withval}" != xyes && test ! -d "${withval}" ; then
      # user may use --with-mpi without an argument, which results in withval
      # being "yes". This case is OK and we simply take no action, as PnetCDF
      # requires MPI compilers and will check them.
      AC_MSG_ERROR(Directory '${withval}' specified in --with-mpi does not exist or is not a directory)
   fi
   MPI_INSTALL=${withval}
   ]
)
if test "x$MPI_INSTALL" = xyes ; then
   # when "--with-mpi" is used
   MPI_INSTALL=
fi

dnl MPI_INSTALL will be referred in UD_MPI_PATH_PROGS and UD_MPI_PATH_PROG
dnl defined in acinclude.m4

AC_ARG_VAR(MPICC,  [MPI C compiler, @<:@default: CC@:>@])
AC_ARG_VAR(MPICXX, [MPI C++ compiler, @<:@default: CXX@:>@])
AC_ARG_VAR(MPIF77, [MPI Fortran 77 compiler, @<:@default: F77@:>@])
AC_ARG_VAR(MPIF90, [MPI Fortran 90 compiler, @<:@default: FC@:>@])
dnl AC_ARG_VAR(CPPFLAGS, [Preprocessor options for C and C++ compilers, e.g. -I<include_dir> if you have headers in a nonstandard directory <include_dir>])
dnl AC_ARG_VAR(CFLAGS, Debugging and optimization options for the C compiler)
dnl AC_ARG_VAR(CXXFLAGS, Debugging and optimization options for the C++ compiler)
dnl AC_ARG_VAR(FFLAGS, Debugging and optimization options for the Fortran 77 compiler)
dnl AC_ARG_VAR(FCFLAGS, Debugging and optimization options for the Fortran 90 compiler)

dnl Check if MPICC, MPICXX, MPIF77, MPIF90 are set by the user.
dnl If not, set MPICC equal to CC. Similarly, for MPICXX, MPIF77, and MPIF90
ac_user_MPICC=$MPICC
ac_user_MPICXX=$MPICXX
ac_user_MPIF77=$MPIF77
ac_user_MPIF90=$MPIF90
if test "x$MPICC"  = x && test "x$CC"  != x ; then ac_user_MPICC=$CC   ; fi
if test "x$MPICXX" = x && test "x$CXX" != x ; then ac_user_MPICXX=$CXX ; fi
if test "x$MPIF77" = x && test "x$F77" != x ; then ac_user_MPIF77=$F77 ; fi
if test "x$MPIF77" = x && test "x$FC"  != x ; then ac_user_MPIF77=$FC  ; fi
if test "x$MPIF90" = x && test "x$F90" != x ; then ac_user_MPIF90=$F90 ; fi
if test "x$MPIF90" = x && test "x$FC"  != x ; then ac_user_MPIF90=$FC  ; fi

CANDIDATE_MPICC="${MPICC} mpicc mpicc_r"
CANDIDATE_MPICXX="${MPICXX} mpicxx mpic++ mpiCC mpcxx mpc++ mpicxx_r mpiCC_r mpcxx_r mpic++_r mpc++_r"
CANDIDATE_MPIF77="${MPIF77} mpif77 mpif77_r mpf77 mpf77_r"
CANDIDATE_MPIF90="${MPIF90} mpif90 mpif90_r mpf90 mpf90_r mpif95 mpif95_r mpf95 mpf95_r"

dnl add GNU MPI compilers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpigcc mpgcc mpigcc_r mpgcc_r"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpig++ mpg++ mpig++_r mpg++_r"
CANDIDATE_MPIF77="$CANDIDATE_MPIF77 mpig77 mpig77_r mpg77 mpg77_r"
CANDIDATE_MPIF90="$CANDIDATE_MPIF90 mpigfortran mpgfortran mpigfortran_r mpgfortran_r"

dnl add IBM MPI compilers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpcc_r mpcc mpixlc_r mpixlc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpCC_r mpCC mpixlcxx_r mpixlcxx mpixlC_r mpixlC"
CANDIDATE_MPIF77="$CANDIDATE_MPIF77 mpixlf77_r mpixlf77"
CANDIDATE_MPIF90="$CANDIDATE_MPIF90 mpixlf90_r mpixlf90"

dnl add IBM BGL MPI compilers
CANDIDATE_MPICC="$CANDIDATE_MPICC blrts_xlc mpxlc_r mpxlc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX blrts_xlC mpxlC_r mpxlC mpixlc++ mpxlcxx mpxlc++ mpxlCC mpixlc++_r mpxlcxx_r mpxlc++_r mpxlCC_r"
CANDIDATE_MPIF77="$CANDIDATE_MPIF77 blrts_xlf mpxlf_r mpxlf mpxlf77 mpxlf77_r mpixlf mpixlf_r"
CANDIDATE_MPIF90="$CANDIDATE_MPIF90 blrts_xlf90 mpxlf90_r mpxlf90 mpxlf95_r mpxlf95 mpixlf95 mpixlf95_r"

dnl add Fujitsu MPI compilers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpifccpx"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpiFCCpx"
CANDIDATE_MPIF77="$CANDIDATE_MPIF77 mpifrtpx"
CANDIDATE_MPIF90="$CANDIDATE_MPIF90 mpifrtpx"

dnl add Cray MPI compiler wrappers
CANDIDATE_MPICC="$CANDIDATE_MPICC cc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX CC"
CANDIDATE_MPIF77="$CANDIDATE_MPIF77 ftn"
CANDIDATE_MPIF90="$CANDIDATE_MPIF90 ftn"

dnl add Intel MPI compiler wrappers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpiicc icc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpiicpc mpiicxx mpiic++ mpiiCC icpc"
CANDIDATE_MPIF77="$CANDIDATE_MPIF77 mpiifort mpiifc mpiif77 ifort"
CANDIDATE_MPIF90="$CANDIDATE_MPIF90 mpiifort mpiifc mpiif95 mpiif90 ifort"

dnl add PGI MPI compiler wrappers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpipgcc mppgcc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpipgCC mppgCC"
CANDIDATE_MPIF77="$CANDIDATE_MPIF77 mpipgf77 mppgf77"
CANDIDATE_MPIF90="$CANDIDATE_MPIF90 mpipgf90 mpipgf95 mppgf90 mppgf95 "

dnl find the full path of MPICC from CANDIDATE_MPICC and MPI_INSTALL
UD_MSG_DEBUG([ac_user_MPICC=$ac_user_MPICC])
if test "x${ac_user_MPICC}" = x ; then
   dnl if MPICC or CC has not been set by users, then search from
   dnl CANDIDATE_MPICC, and find the full path of MPICC
   UD_MPI_PATH_PROGS([MPICC], [$CANDIDATE_MPICC])
else
   dnl check whether user specified MPICC is valid
   UD_MPI_PATH_PROG([MPICC], [$ac_user_MPICC])
fi
UD_MSG_DEBUG([MPICC=$MPICC])

if test "x${MPICC}" = x ; then
   if test "x$ac_user_MPICC" = x ; then
      ERR_MSG="No MPI C compiler can be found"
   else
      ERR_MSG="Specified MPI C compiler \"$ac_user_MPICC\" cannot be found"
   fi
   if test "x$MPI_INSTALL" != x ; then
      ERR_MSG="$ERR_MSG under $MPI_INSTALL"
   fi
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
     $ERR_MSG
     PnetCDF requires a working MPI C compiler. Please specify the location
     of an MPI C compiler, either in the MPICC environment variable
     (not CC variable) or through --with-mpi configure flag. Abort.
   -----------------------------------------------------------------------])
fi
UD_MSG_DEBUG([MPICC=$MPICC])

dnl Now MPICC is configured, do some basic compiler tests
CC=${MPICC}
AC_PROG_CC
dnl AM_PROG_CC_C_O

dnl get compiler vendor in ax_cv_c_compiler_vendor (e.g. gnu, intel)
AX_COMPILER_VENDOR

if test "x${ax_cv_c_compiler_vendor}" = xgnu ; then
   dnl gcc command-line option "-dumpversion" can also show version.
   AC_MSG_CHECKING([whether gcc version is greater than 10.0.0])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#if __GNUC__ < 10
#error gcc version < 10.0.0
#endif
   ]])], [gcc_ge_10=yes], [gcc_ge_10=no])
   AC_MSG_RESULT([$gcc_ge_10])
fi

dnl get base compiler command of MPI C compiler wrapper in
dnl ac_cv_mpi_compiler_base_MPICC (e.g. /usr/bin/gcc)
MPI_COMPILER_BASE(MPICC)

dnl Set output variable CPP to a command that runs the C preprocessor.
dnl Some C compilers require -E to be used as C preprocessor.
AC_PROG_CPP

dnl this call needs at least autoconf version 2.60
dnl AC_USE_SYSTEM_EXTENSIONS

dnl enable large file support
AC_SYS_LARGEFILE

dnl UD_PROG_CC_MAKEDEPEND

dnl AC_HEADER_STDC
dnl AC_CHECK_HEADERS([malloc.h])

dnl check if MPICC works for basic MPI call: MPI_Comm_rank()
AC_CHECK_FUNC([MPI_Comm_rank], [],
   dnl maybe -lmpi is needed at link stage
   [AC_SEARCH_LIBS([MPI_Comm_rank], [mpi mpich], [],
                 [AC_MSG_ERROR([
   -----------------------------------------------------------------------
     Invalid MPI compiler specified or detected: "${MPICC}"
     A working MPI C compiler is required. Please specify the location
     of one either in the MPICC environment variable (not CC variable) or
     through --with-mpi configure flag. Abort.
   -----------------------------------------------------------------------])
])])

AC_CHECK_FUNC([MPI_File_open], [],
   dnl maybe -lmpio is needed at link stage
   [AC_SEARCH_LIBS([MPI_File_open], [mpi mpio], [],
                   [AC_MSG_ERROR([
   -----------------------------------------------------------------------
     The underneath MPI implementation does not support MPI-IO.
     PnetCDF requires MPI-IO support to work properly. Abort.
   -----------------------------------------------------------------------])
])])

dnl check the base of MPICC compiler
AM_CONDITIONAL(MPICC_IS_XLC,  [test "x${ax_cv_c_compiler_vendor}" = xibm])
AM_CONDITIONAL(MPICC_IS_PGCC, [test "x${ax_cv_c_compiler_vendor}" = xportland])
AM_CONDITIONAL(MPICC_IS_FCCPX,[test "x${ax_cv_c_compiler_vendor}" = xfujitsu])

dnl print out MPI version and vendor information
CHECK_MPI_VERSION

dnl PnetCDF requires MPI-IO support and MPI-IO was first introduced in MPI 2.0
if test "$mpi_version" -lt "2" ; then
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
     PnetCDF requires the MPI library that supports MPI standard version
     2.0 and later. The supplied MPI library only supports MPI $mpi_version.
     Abort.
   -----------------------------------------------------------------------])
fi

dnl check whether MPICC is built using a C++ compiler
AC_MSG_CHECKING([whether $MPICC is a wrapper of a C++ compiler])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [int class=0;])],
                  [AC_MSG_RESULT(no)],
                  [AC_MSG_ERROR([
   -----------------------------------------------------------------------
     PnetCDF requires the supplied MPI C compiler is not a wrapper of a
     C++ compiler. MPI C compiler is $MPICC
     Abort.
   -----------------------------------------------------------------------])]
)

dnl compute canonical system types
AC_CANONICAL_BUILD
AC_CANONICAL_HOST

dnl AM_PROG_AR is first recognized/traced in autoconf 2.68
dnl Use it only if we want support for unusual archivers such as Microsoft
dnl lib.. AM_PROG_AR must be called before LT_INIT
AM_PROG_AR
dnl UD_PROG_AR()
dnl UD_PROG_NM()
dnl We could use the PAC check for ranlib (it also makes sure that ranlib
dnl works, which is not always true, particularly when GNU tools are
dnl installed on a system that does not have (or need) ranlib
dnl libtoolize: `AC_PROG_RANLIB' is rendered obsolete by `LT_INIT'
dnl AC_PROG_RANLIB

dnl libtool v2.4.6 was released in 2015-02-15
dnl Travis CI only has v2.4.2
LT_PREREQ([2.4.2])
dnl LT_INIT([dlopen disable-shared])
dnl LT_INIT([dlopen])
LT_INIT([disable-shared])

dnl check MPI C++ compiler
AC_ARG_ENABLE(cxx,
    [AS_HELP_STRING([--disable-cxx],
                    [Turn off support for the C++ interface,
                     if you only need the C interface. @<:@default: enabled@:>@])],
    [enable_cxx=${enableval}], [enable_cxx=auto]
)
UD_MSG_DEBUG(enable_cxx=$enable_cxx)

dnl If MPICXX or CXX has not been set by users, then
dnl search from CANDIDATE_MPICXX and find the full path of MPICXX
if test "x${enable_cxx}" != xno ; then
   if test "x${ac_user_MPICXX}" = x ; then
      dnl if MPICXX or CXX has not been set by users, then search from
      dnl CANDIDATE_MPICXX, and find the full path of MPICXX
      UD_MPI_PATH_PROGS([MPICXX], [$CANDIDATE_MPICXX])
   else
      dnl check whether user specified MPICXX is valid
      UD_MPI_PATH_PROG([MPICXX], [$ac_user_MPICXX])
   fi
   UD_MSG_DEBUG([MPICXX=$MPICXX])

   if test "x${MPICXX}" = x ; then
      if test "x$ac_user_MPICXX" = x ; then
         ERR_MSG="No MPI C++ compiler can be found"
      else
         ERR_MSG="Specified MPI C++ compiler \"$ac_user_MPICXX\" cannot be found"
      fi
      if test "x$MPI_INSTALL" != x ; then
         ERR_MSG="$ERR_MSG under $MPI_INSTALL"
      fi
      if test "x${enable_cxx}" = xyes ; then
         dnl --enable-cxx is explicitly set at command line
         AC_MSG_ERROR([
   -----------------------------------------------------------------------
     $ERR_MSG.
     Please specify the location of the MPI C++ compiler, either in the
     MPICXX environment variable or the --with-mpi configure flag. Abort.
   -----------------------------------------------------------------------])
      else
         dnl enable_cxx is auto, i.e. --enable-cxx is not set at command line
         AC_MSG_WARN([
   -----------------------------------------------------------------------
     $ERR_MSG.
     Thus the C++ feature is disabled.
     The location of the MPI C++ compiler can be specified either in the
     MPICXX environment variable or through the --with-mpi configure flag.
   -----------------------------------------------------------------------])
      fi
   fi
fi

if test "x${enable_cxx}" = xno || test "x${MPICXX}" = x ; then
   dnl when cxx is explicitly disabled or no MPICXX can be found
   has_mpicxx=no
else
   has_mpicxx=yes
   CXX=${MPICXX}
fi

dnl Must invoke AC_PROG_CXX here, because AC_PROG_CXX cannot be called
dnl conditionally in Automake, otherwise got "error: conditional
dnl "am__fastdepCXX" was never defined." So, if has_mpicxx is no, then
dnl AC_PROG_CXX will search for g++. However, g++ will not be used at all.
AC_PROG_CXX

dnl autoconf 2.59 has not yet implemented AC_PROG_CXX_C_O
dnl AC_PROG_CXX_C_O

if test "x${has_mpicxx}" = xyes ; then
   dnl test if MPICXX can compile an MPI-IO program
   AC_LANG_PUSH(C++)
   dnl get compiler vendor in ax_cv_cxx_compiler_vendor (e.g. gnu, intel)
   AX_COMPILER_VENDOR
   dnl get base compiler command of MPI C++ compiler wrapper in
   dnl ac_cv_mpi_compiler_base_MPICXX (e.g. /usr/bin/g++)
   MPI_COMPILER_BASE(MPICXX)
   AC_CHECK_FUNC([MPI_File_close], [],
      dnl maybe -lmpi++ is needed at link stage
      [AC_SEARCH_LIBS([MPI_File_close], [mpi mpio mpi++ mpichcxx mpi_cxx], [],
                      [has_mpicxx=no])])
   AC_LANG_POP(C++)
   if test "x$has_mpicxx" = xno ; then
      if test "x${enable_cxx}" = xyes ; then
         AC_MSG_ERROR([
   -----------------------------------------------------------------------
     Explicitly requesting C++ feature, but "${MPICXX}"
     is not a working MPI C++ compiler. Abort.
   -----------------------------------------------------------------------])
      else dnl in case enable_cxx is auto
         AC_MSG_WARN([
   -----------------------------------------------------------------------
     "${MPICXX}" is not a working MPI C++ compiler.
     Thus, the C++ feature is disabled.
   -----------------------------------------------------------------------])
      fi
   else
      UD_CHECK_MPI_CPP_SEEK_SET
      UD_MSG_DEBUG(ac_cv_CHECK_MPI_CPP_SEEK_SET=$ac_cv_CHECK_MPI_CPP_SEEK_SET)
      dnl When using older version Intel compilers 4.x, SEEK_SET will be
      dnl reported as redefined. We need to add the following C++ preprocessor
      dnl flags. See doc/README.INTEL
      if test "x${ac_cv_CHECK_MPI_CPP_SEEK_SET}" = xyes ; then
         SEEK_SET_REDEFINED=yes
      fi
   fi
fi
if test "x${has_mpicxx}" = xyes ; then
   ENABLE_CXX=1
   AC_DEFINE(ENABLE_CXX)
else
   ENABLE_CXX=0
fi
AC_SUBST(ENABLE_CXX)

AC_SUBST(has_mpicxx)dnl for src/utils/pnetcdf-config.in
UD_MSG_DEBUG(has_mpicxx=$has_mpicxx)
AM_CONDITIONAL(HAS_MPICXX, [test x$has_mpicxx = xyes])
AM_CONDITIONAL(SEEK_SET_REDEFINED, [test x$ac_cv_CHECK_MPI_CPP_SEEK_SET = xyes])

dnl Note this must be done after the type of C compiler is determined
AC_ARG_ENABLE(strict,
    [AS_HELP_STRING([--enable-strict],
                    [Turn on strict debugging with gcc. @<:@default: disabled@:>@])],
    [enable_strict=${enableval}], [enable_strict=no]
)

if test "x${enable_strict}" = xyes; then
   if test "x${GCC}" = xyes; then
      CFLAGS="$CFLAGS -Wall -Wstrict-prototypes -Wmissing-prototypes -Wundef -Wpointer-arith -Wbad-function-cast"
   else
      AC_MSG_WARN([--enable-strict is for GNU compiler only])
      enable_strict=no
   fi
fi

AC_ARG_ENABLE(fortran,
    [AS_HELP_STRING([--disable-fortran],
                    [Turn off support for the Fortran interface,
                     if you only need the C interface. @<:@default: enabled@:>@])],
    [enable_fortran=${enableval}], [enable_fortran=auto]
)

if test "x${enable_fortran}" = xno ; then
   has_fortran=no
else
   has_fortran=yes
   if test "x${ac_user_MPIF77}" = x ; then
      dnl if MPIF77 or F77 has not been set by users, then search from
      dnl CANDIDATE_MPIF77, and find the full path of MPIF77
      UD_MPI_PATH_PROGS([MPIF77], [$CANDIDATE_MPIF77])
   else
      dnl check whether user specified MPIF77 is valid
      UD_MPI_PATH_PROG([MPIF77], [$ac_user_MPIF77])
   fi

   if test "x${MPIF77}" = x ; then
      has_fortran=no
      if test "x$ac_user_MPIF77" = x ; then
         ERR_MSG="No MPI Fortran 77 compiler can be found"
      else
         ERR_MSG="Specified MPI Fortran 77 compiler \"$ac_user_MPIF77\" cannot be found"
      fi
      if test "x$MPI_INSTALL" != x ; then
         ERR_MSG="$ERR_MSG under $MPI_INSTALL"
      fi
      if test "x${enable_fortran}" = xyes ; then
         dnl --enable-fortran is explicitly set at command line
         AC_MSG_ERROR([
   -----------------------------------------------------------------------
     $ERR_MSG.
     Please specify the location of the MPI Fortran 77 compiler, either in the
     MPIF77 environment variable or the --with-mpi configure flag. Abort.
   -----------------------------------------------------------------------])
      else
         dnl enable_fortran is auto, i.e. --enable-fortran is not set at command line
         AC_MSG_WARN([
   -----------------------------------------------------------------------
     $ERR_MSG.
     Thus, the Fortran feature is disabled.
     The location of the MPI Fortran 77 compiler can be specified either in the
     MPIF77 environment variable or through the --with-mpi configure flag.
   -----------------------------------------------------------------------])
      fi
   fi
fi

if test "x${has_fortran}" = xyes ; then
   if test "x${ac_user_MPIF90}" = x ; then
      dnl if MPIF90 or F90 has not been set by users, then search from
      dnl CANDIDATE_MPIF90, and find the full path of MPIF90
      UD_MPI_PATH_PROGS([MPIF90], [$CANDIDATE_MPIF90])
   else
      dnl check whether user specified MPIF90 is valid
      UD_MPI_PATH_PROG([MPIF90], [$ac_user_MPIF90])
   fi

   if test "x${MPIF90}" = x ; then
      has_fortran=no
      if test "x$ac_user_MPIF90" = x ; then
         ERR_MSG="No MPI Fortran 90 compiler can be found"
      else
         ERR_MSG="Specified MPI Fortran 90 compiler \"$ac_user_MPIF90\" cannot be found"
      fi
      if test "x$MPI_INSTALL" != x ; then
         ERR_MSG="$ERR_MSG under $MPI_INSTALL"
      fi
      if test "x${enable_fortran}" = xyes ; then
         dnl --enable-fortran is explicitly set at command line
         AC_MSG_ERROR([
   -----------------------------------------------------------------------
     $ERR_MSG.
     Please specify the location of the MPI Fortran 90 compiler, either in the
     MPIF90 environment variable or the --with-mpi configure flag. Abort.
   -----------------------------------------------------------------------])
      else
         dnl enable_fortran is auto, i.e. --enable-fortran is not set at command line
         AC_MSG_WARN([
   -----------------------------------------------------------------------
     $ERR_MSG.
     Thus, the Fortran feature is disabled.
     The location of the MPI Fortran 90 compiler can be specified either in the
     MPIF90 environment variable or through the --with-mpi configure flag.
   -----------------------------------------------------------------------])
      fi
   fi
fi
UD_MSG_DEBUG(has_fortran=${has_fortran})
UD_MSG_DEBUG([MPIF77=$MPIF77])
UD_MSG_DEBUG([MPIF90=$MPIF90])

if test "x${has_fortran}" = xyes ; then

   dnl Check if MPIF77 is a valid MPI compiler
   F77=${MPIF77}
   AC_PROG_F77
   dnl FFLAGS is set in AC_PROG_F77

   AC_F77_LIBRARY_LDFLAGS

   AC_MSG_CHECKING([whether $MPIF77 is a valid MPI compiler])
   AC_LANG_PUSH([Fortran 77])
   dnl check if can use mpif.h
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
           include "mpif.h"
           integer err, rank
           call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)]])],
      [valid_mpif77=yes],[valid_mpif77=no])
   AC_LANG_POP([Fortran 77])
   AC_MSG_RESULT($valid_mpif77)

   if test "x${valid_mpif77}" = xno ; then
      has_fortran=no
      if test "x${enable_fortran}" = xyes ; then
         AC_MSG_ERROR([
   -----------------------------------------------------------------------
     Invalid MPI Fortran 77 compiler: "${MPIF77}"
     A working MPI Fortran 77 compiler is required. Please specify the
     location of a valid MPI Fortran 77 compiler, either in the MPIF77
     environment variable or through --with-mpi configure flag. Abort.
   -----------------------------------------------------------------------])
      else
         AC_MSG_WARN([
   -----------------------------------------------------------------------
     "${MPIF77}" is not a working MPI Fortran 77 compiler.
     Thus, the Fortran feature is disabled.
   -----------------------------------------------------------------------])
      fi
   fi
fi

mpi_mod=no
if test "x${has_fortran}" = xyes ; then

   dnl Check if MPIF90 is a valid MPI compiler
   FC=${MPIF90}
   AC_PROG_FC
   dnl FCFLAGS is set in AC_PROG_FC

   dnl FCLIBS_save="$FCLIBS"
   dnl FCLIBS=""
   AC_FC_LIBRARY_LDFLAGS
   dnl UD_MSG_DEBUG([before FCLIBS=$FCLIBS])
   dnl The autoconf macro for finding FCLIBS sometimes makes mistakes
   dnl (particularly with the Fujitsu frt compiler).  This next step
   dnl first sees if the FCLIBS is valid with the Fortran compiler
   dnl This also happens to Solaris Studio Fortran compilers
   dnl AC_PROG_FC_FCLIBS_VALID
   dnl Now see if FCLIBS works with the C compiler (remove invalid ones)
   dnl PAC_PROG_FC_CHECK_FCLIBS

   dnl replace FLIBS and F90LIBS with FCLIBS
   dnl UD_MSG_DEBUG([after FCLIBS=$FCLIBS])
   dnl FLIBS="$FCLIBS"
   dnl F90LIBS="$FCLIBS"

   AC_MSG_CHECKING([whether mpi.mod is available])
   AC_LANG_PUSH([Fortran])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[        use mpi]])],
                     [mpi_mod=yes], [mpi_mod=no]
   )
   AC_MSG_RESULT($mpi_mod)
   AC_MSG_CHECKING([whether $MPIF90 is a valid MPI compiler])
   if test "x${mpi_mod}" = xyes ; then
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
           use mpi
           integer err, rank
           call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)]])],
         [valid_mpif90=yes],[valid_mpif90=no]
      )
   else
      AC_MSG_WARN([Fortran module mpi.mod is not available or invalid. Use mpif.h instead.])
      dnl mpi.mod is not available, check if can use mpif.h
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
           include "mpif.h"
           integer err, rank
           call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)]])],
         [valid_mpif90=yes],[valid_mpif90=no]
      )
   fi
   AC_LANG_POP([Fortran])
   AC_MSG_RESULT($valid_mpif90)

   if test "x${valid_mpif90}" = xno ; then
      has_fortran=no
      if test "x${enable_fortran}" = xyes ; then
         AC_MSG_ERROR([
   -----------------------------------------------------------------------
     Invalid MPI Fortran 90 compiler: "${MPIF90}"
     A working MPI Fortran 90 compiler is required. Please specify the
     location of a valid MPI Fortran 90 compiler, either in the MPIF90
     environment variable or through --with-mpi configure flag. Abort.
   -----------------------------------------------------------------------])
      fi
      AC_MSG_WARN([
   -----------------------------------------------------------------------
     "${MPIF90}" is not a working MPI Fortran 90 compiler.
     Thus, the Fortran feature is disabled.
   -----------------------------------------------------------------------])
   else
      if test "x${mpi_mod}" = xyes ; then
         USE_MPIF_HEADER="use mpi, only: MPI_OFFSET_KIND"
      else
         USE_MPIF_HEADER="include \"mpif.h\""
      fi
      AC_SUBST(USE_MPIF_HEADER)
   fi
fi
if test "x${has_fortran}" = xyes ; then
   ENABLE_FORTRAN=1
   AC_DEFINE(ENABLE_FORTRAN)
else
   ENABLE_FORTRAN=0
fi
AC_SUBST(ENABLE_FORTRAN)

AC_SUBST(has_fortran)dnl for src/utils/pnetcdf-config.in
AM_CONDITIONAL(HAS_FORTRAN, [test x$has_fortran = xyes])
AM_CONDITIONAL(HAVE_MPI_MOD, [test x$mpi_mod = xyes])

if test "x${has_fortran}" = xyes ; then
   if test "x${enable_strict}" = xyes && test "x${GFC}" = xyes; then
      dnl not all Fortran compilers recognize -Wall
      FFLAGS="$FFLAGS -Wall"
      FCFLAGS="$FCFLAGS -Wall"
   fi
   dnl AC_SUBST(FLIBS)
   dnl AC_SUBST(FCLIBS)
   dnl AC_SUBST(F90LIBS)
   dnl AC_SUBST(FLDFLAGS)
   dnl AC_SUBST(F90LDFLAGS)
fi

ac_f77_support_freeform=no
if test "x${has_fortran}" = xyes ; then
   FC_saved=${FC}
   FC=${MPIF77}
   dnl customized AC_FC_FREEFORM: to just get ac_cv_fc_freeform without
   dnl appending ac_cv_fc_freeform to FCFLAGS
   dnl ac_cv_fc_freeform is the flag for enabling Fortran free form
   dnl we use this flag only in testing programs
   UD_FC_FREEFORM
   if test "x${ac_cv_fc_freeform}" != xunknown ; then
      ac_f77_support_freeform=yes
      FFREEFORMFLAG=${ac_cv_fc_freeform}
   fi
   UD_MSG_DEBUG([FFREEFORMFLAG=$FFREEFORMFLAG])
   AC_SUBST(FFREEFORMFLAG)
   FC=${FC_saved}
fi
AM_CONDITIONAL(HAVE_F77_SUPPORT_FREEFORM, [test x$ac_f77_support_freeform = xyes])

if test "x${has_fortran}" = xyes ; then
   dnl GNU Fortran compiler automatically invokes preprocessor for files with
   dnl extension .F and .F90. To manually invoke the preprocessor on any file,
   dnl use compiler flag -cpp. To disable, use -nocpp.
   dnl For other compilers, we need to find the compile flags for Fortran
   dnl preprocessor.
   dnl AC_FC_PP_SRCEXT is first introduced in autoconf 2.69
   dnl steal AC_FC_PP_SRCEXT from autoconf 2.69 to make UD_FC_PP_SRCEXT
   dnl UD_FC_PP_SRCEXT([F])    dnl sets ac_cv_fc_pp_srcext_F
   dnl UD_FC_PP_SRCEXT([F90])  dnl sets ac_cv_fc_pp_srcext_F90
   AC_FC_PP_SRCEXT([F])    dnl sets ac_cv_fc_pp_srcext_F
   AC_FC_PP_SRCEXT([F90])  dnl sets ac_cv_fc_pp_srcext_F90

   dnl compiler command-line define preprocessor flag, result in ${FC_DEFINE}
   dnl Not all Fortran compilers use -D
   dnl UD_FC_PP_DEFINE
   dnl AC_FC_PP_DEFINE is first introduced in autoconf 2.69
   AC_FC_PP_DEFINE

   dnl check compiler flags for file extensions in .f .F .f90 .F90
   AC_FC_SRCEXT([f])
   AC_FC_SRCEXT([F])
   AC_FC_SRCEXT([f90])
   AC_FC_SRCEXT([F90])

   AC_LANG_PUSH([Fortran])
   dnl get compiler vendor in ax_cv_fc_compiler_vendor (e.g. gnu, intel)
   AX_COMPILER_VENDOR
   AC_LANG_POP([Fortran])

   dnl get base compiler command of MPI Fortran compiler wrapper in
   dnl ac_cv_mpi_compiler_base_MPIF90 (e.g. /usr/bin/gfortran)
   MPI_COMPILER_BASE(MPIF90)

   dnl customized AC_FC_FIXEDFORM: to just get ac_cv_fc_fixedform without
   dnl appending ac_cv_fc_fixedform to FCFLAGS
   dnl ac_cv_fc_fixedform flag is for enabling Fortran fixed form.
   dnl we use this flag only in testing programs
   UD_FC_FIXEDFORM
   FFIXEDFORMFLAG=${ac_cv_fc_fixedform}
   UD_MSG_DEBUG([FC=$FC FFIXEDFORMFLAG=$FFIXEDFORMFLAG])
   AC_SUBST(FFIXEDFORMFLAG)

   dnl Checking for Fortran types also determines the Fortran name mangling
   dnl and places the value into FCALLSCSUB as the C name corresponding
   dnl to the Fortran name SUB
   AC_FC_FUNC(sub, [FCALLSCSUB])

   dnl determine the correct name mapping
   case $FCALLSCSUB in
       SUB)
           AC_DEFINE(F77_NAME_UPPER,,[Define if Fortran names are uppercase])
           ;;
       sub_)
           dnl This is the hard case.  gcc uses one _ unless the name includes
           dnl an underscore, in which case it gets two trailing underscores.
           dnl Use essentially the same configure code that the original configure
           dnl used to determine SUB
           AC_MSG_CHECKING([for C-equivalent to Fortran routine "SUB_A"])
                   dnl "
           AC_FC_FUNC(sub_a, [FCALLSCSUBA])
           AC_MSG_RESULT($FCALLSCSUBA)
           case $FCALLSCSUBA in
               sub_a__)
                   AC_DEFINE(F77_NAME_LOWER_2USCORE,,[Define if Fortran names are lower case with two trailing underscore2])
                   ;;
               sub_a_)
                   AC_DEFINE(F77_NAME_LOWER_USCORE,,[Define if Fortran names are lower case with one trailing underscore])
                   ;;
               *)
                   AC_MSG_WARN([Unrecognized Fortran name mapping])
                   ;;
           esac
           ;;
       sub)
           AC_DEFINE(F77_NAME_LOWER,,[Define if Fortran names are lower case])
           ;;
       *)
           AC_MSG_WARN([Unrecognized Fortran name mapping])
           ;;
   esac

   dnl Some Fortran 77 compilers, such as pgf77, do not allow "_8" modifier,
   dnl because _8 modifier is a Fortran 90 feature
   dnl UD_CHECK_F77_IS_PGF77
   dnl UD_MSG_DEBUG([ac_cv_mpif77_is_PGF77=$ac_cv_mpif77_is_PGF77])

   dnl check Fortran parameter modifier for 8-byte integer type
   dnl We need this to set the max constants for UINT, INT64, and UINT64
   UD_FC_CONSTANT_MODIFIER
   UD_MSG_DEBUG([ac_cv_fc_constant_modifier=$ac_cv_fc_constant_modifier])
   PNF_INT8_MODIFIER=""
   if test "x${ac_cv_fc_constant_modifier}" = xnone ; then
      PNF_FILL_UINT=4294967295
      PNF_FILL_INT64=-9223372036854775806
      PNF_FILL_UINT64=18446744073709551614
      PNF_X_UINT_MAX=4294967295
      PNF_X_INT8_MIN=-9223372036854775807
      PNF_X_INT8_MAX=9223372036854775807
      PNF_X_UINT8_MAX=18446744073709551615
   else
      if test "x${ac_cv_fc_constant_modifier}" = xEightByteInt ; then
         PNF_INT8_MODIFIER="      integer, parameter :: EightByteInt = selected_int_kind(18)"
      fi
      PNF_FILL_UINT=4294967295_${ac_cv_fc_constant_modifier}
      PNF_FILL_INT64=-9223372036854775806_${ac_cv_fc_constant_modifier}
      PNF_FILL_UINT64=18446744073709551614_${ac_cv_fc_constant_modifier}
      PNF_X_UINT_MAX=4294967295_${ac_cv_fc_constant_modifier}
      PNF_X_INT8_MIN=-9223372036854775807_${ac_cv_fc_constant_modifier}
      PNF_X_INT8_MAX=9223372036854775807_${ac_cv_fc_constant_modifier}
      PNF_X_UINT8_MAX=18446744073709551615_${ac_cv_fc_constant_modifier}
   fi
   AC_SUBST(PNF_INT8_MODIFIER)
   AC_SUBST(PNF_FILL_UINT)
   AC_SUBST(PNF_FILL_INT64)
   AC_SUBST(PNF_FILL_UINT64)
   AC_SUBST(PNF_X_UINT_MAX)
   AC_SUBST(PNF_X_INT8_MIN)
   AC_SUBST(PNF_X_INT8_MAX)
   AC_SUBST(PNF_X_UINT8_MAX)
fi

ac_cv_prog_f90_uppercase_mod=no
if test "x${has_fortran}" = xyes ; then
   dnl
   dnl Check Fortran module file extension
   dnl
   dnl UD_FC_MODULE_EXTENSION
   AC_FC_MODULE_EXTENSION
   dnl UD_FC_MODULE_EXTENSION defines FC_MODEXT
   if test "x${FC_MODEXT}" = x ; then
      AC_MSG_ERROR([cannot determine Fortran module file extension!])
   fi
   UD_MSG_DEBUG([FC_MODEXT=$FC_MODEXT])

   dnl AC_FC_MODULE_FLAG was first introduced in autoconf 2.69
   dnl AC_FC_MODULE_FLAG/UD_FC_MODULE_FLAG defines FC_MODINC
   dnl However, AC_FC_MODULE_FLAG has a bug that mistakenly sets FC_MODINC to
   dnl -M when Fujitsu frtpx is used. UD_FC_MODULE_FLAG fixes this problem.
   UD_FC_MODULE_FLAG
   UD_MSG_DEBUG([FC_MODINC=$FC_MODINC])

   dnl AC_FC_MODULE_OUTPUT_FLAG
   dnl UD_FC_MODULE_OUTPUT_FLAG
   dnl UD_FC_MODULE_OUTPUT_FLAG defines FC_MODOUT
   dnl UD_MSG_DEBUG([FC_MODOUT=$FC_MODOUT])

   dnl
   dnl Below is to check if a Fortran compiler produces module files with upper
   dnl case file name, e.g. PNETCDF.mod. However, this does not work for Mac
   dnl OSX file system which is case insensitive
   dnl
   UD_PROG_FC_UPPERCASE_MOD
fi
AM_CONDITIONAL(UPPER_CASE_MOD, [test x$ac_cv_prog_f90_uppercase_mod = xyes])

if test "x${has_fortran}" = xyes && test "x${ax_cv_fc_compiler_vendor}" = xnag ; then
   # Special treatment for NAG Fortran compiler
   # Command 'nagfor -dryrun dummy.f90' can reveal the C linker nagfor uses.
   if test "x${ac_cv_mpi_compiler_base_MPICC}" != x ; then
      UD_MSG_DEBUG([ac_cv_mpi_compiler_base_MPICC=$ac_cv_mpi_compiler_base_MPICC])
      str_found=`echo "${FFLAGS}" | ${EGREP} -- "-Wc="`
      if test "x$?" != x0 ; then
         AC_MSG_WARN([Add to FFLAGS -Wc=$ac_cv_mpi_compiler_base_MPICC])
         AS_VAR_APPEND([FFLAGS], [" -Wc=$ac_cv_mpi_compiler_base_MPICC"])
      fi
      str_found=`echo "${FCFLAGS}" | ${EGREP} -- "-Wc="`
      if test "x$?" != x0 ; then
         AC_MSG_WARN([Add to FCFLAGS -Wc=$ac_cv_mpi_compiler_base_MPICC])
         AS_VAR_APPEND([FCFLAGS], [" -Wc=$ac_cv_mpi_compiler_base_MPICC"])
      fi
      unset str_found
   fi

   # Add -mismatch to downgrade consistency checking of procedure argument
   # lists so that mismatches produce warning messages instead of error
   # messages
   ACX_F77_MISMATCH([AS_VAR_APPEND([FFLAGS], [" $acx_cv_f77_mismatch_flag"])])
   ACX_FC_MISMATCH([AS_VAR_APPEND([FCFLAGS], [" $acx_cv_fc_mismatch_flag"])])
fi

if test "x$enable_shared" = xyes ; then
   dnl Call LT_OUTPUT to produce file ./libtool which is used later in
   dnl. LT_AC_CHECK_SHLIB. Note LT_OUTPUT must be called after all compilers
   dnl checked, such as AC_PROG_CC, AC_PROG_CXX, AC_PROG_FC
   LT_OUTPUT

   dnl Check whether MPICC supports shared libraries.
   LT_MPI_CHECK_SHLIB
   if test "x$ac_cv_lt_mpi_check_shlib" = xno ; then
      AC_MSG_ERROR([
   -----------------------------------------------------------------------
     configure option --enable-shared is set, but the MPI library is not
     built with shared library support. Either remove option --enable-shared
     or use an MPI library built with shared library support. Abort.
   -----------------------------------------------------------------------])
   fi
fi

dnl AC_PROG_INSTALL

dnl have_yacc_lex=no
dnl Starting from PnetCDF 1.5.0, yacc and lex is no longer needed
dnl see comments in src/utils/ncmpigen/Makefile.in for build rules for
dnl ncmpigenyy.c and ncmpigentab.c. If rebuild is desired, uncomment
dnl below checking to check availability of yacc/lex/bison
dnl AC_PROG_YACC
dnl dnl if neither bison nor byacc is found, YACC will be set to yacc
dnl have_yacc_lex=yes
dnl if test "x$YACC" = xyacc; then
dnl    AC_CHECK_PROGS(YACC_PATH, yacc)
dnl    if test "x$YACC_PATH" = x; then
dnl       dnl cannot find bison or yacc required to build ncmpigentab.c
dnl       have_yacc_lex=no
dnl    fi
dnl fi
dnl
dnl AC_PROG_LEX
dnl if (test "x$LEX" != xflex) && (test "x$LEX" != xlex) ; then
dnl    dnl cannot find flex or lex required to build ncmpigenyy.c
dnl    have_yacc_lex=no
dnl fi
dnl AM_CONDITIONAL(HAVE_YACC_LEX, [test x$have_yacc_lex = xyes])

dnl AC_PROG_LN_S
dnl AC_PROG_MAKE_SET
dnl AC_SUBST(SET_MAKE)

dnl YACC and LEX are required to build PnetCDF utility tool ncmpigen
dnl if configure finds bison then YACC is set to bison -y, so we need to clean
dnl up the output a bit before testing
dnl Below checks commands yacc and lex availability under PATH. However, this
dnl checking is redundant, as AC_PROG_* did that already
dnl YACC_CMD="${YACC% *}"
dnl AC_CHECK_PROG([yacc_cmd], [${YACC_CMD}], [yes], [no])
dnl AC_CHECK_PROG([lex_cmd],  [${LEX}],      [yes], [no])
dnl if test "x${yacc_cmd}" = no ; then
dnl     AC_MSG_ERROR([could not find bison/yacc required by PnetCDF])
dnl fi
dnl if test "x${lex_cmd}" = no ; then
dnl     AC_MSG_ERROR([could not find flex/lex required by PnetCDF])
dnl fi

UD_PROG_M4
M4FLAGS="$M4FLAGS -DPNETCDF"
M4FFLAGS="$M4FFLAGS -DPNETCDF"

dnl Check for <stdbool.h> that conforms to C99 requirements
dnl this is also for using bool type in utf8proc.h/utf8proc.c to support
dnl special characters in CDF-2 and CDF-5
AC_HEADER_STDBOOL

dnl AC_C_CONST
AC_C_INLINE

dnl we do not use struct stat yet
dnl AC_CHECK_MEMBERS([struct stat.st_blksize])

UD_CHECK_IEEE

dnl cross compile fails with undefined reference to rpl_realloc and rpl_malloc
dnl AC_FUNC_MALLOC
dnl AC_FUNC_REALLOC

dnl below checks availability of a bunch C functions, but we have not yet
dnl implemented alternative calls
dnl AC_FUNC_ERROR_AT_LINE
dnl AC_FUNC_MEMCMP
dnl AC_FUNC_STRTOD
dnl AC_FUNC_VPRINTF
dnl AC_CHECK_FUNCS([memset setlocale sqrt strchr strrchr strtol])
dnl AC_CHECK_LIB([m], [tanh])
dnl UD_CHECK_LIB_MATH

dnl When using gcc based compiler with -ansi flag, AC_CHECK_FUNCS can still
dnl find strdup, but AC_CHECK_DECL cannot. So we check with AC_CHECK_DECL
dnl first and then check AC_CHECK_FUNCS.
dnl AC_CHECK_FUNCS([strdup strerror access unlink])
AC_CHECK_DECL([strdup], [], [], [[#include <string.h>]])
if test "x$ac_cv_have_decl_strdup" = xyes ; then
   AC_CHECK_FUNCS([strdup])
else
   # strdup is not declared, we add our own implementation of strdup()
   # do not call AC_REPLACE_FUNCS, as it calls AC_CHECK_FUNCS
   AC_CONFIG_LIBOBJ_DIR([src/libs])
   AC_LIBSOURCE([strdup.c])
   AC_LIBOBJ([strdup])
fi
AC_CHECK_DECL([strerror], [], [], [[#include <string.h>]])
if test "x$ac_cv_have_decl_strerror" = xyes ; then
   AC_CHECK_FUNCS([strerror])
fi
AC_CHECK_DECL([access], [], [], [[#include <unistd.h>]])
if test "x$ac_cv_have_decl_access" = xyes ; then
   AC_CHECK_FUNCS([access])
fi
AC_CHECK_DECL([truncate], [], [], [[#include <unistd.h>
#include <sys/types.h>]])
if test "x$ac_cv_have_decl_truncate" = xyes ; then
   AC_CHECK_FUNCS([truncate])
fi
AC_CHECK_DECL([unlink], [], [], [[#include <unistd.h>]])
if test "x$ac_cv_have_decl_unlink" = xyes ; then
   AC_CHECK_FUNCS([unlink])
fi
AC_CHECK_DECL([strcasecmp], [], [], [[#include <strings.h>]])
if test "x$ac_cv_have_decl_strcasecmp" = xyes ; then
   AC_CHECK_FUNCS([strcasecmp])
else
   # strcasecmp is not declared, we add our own implementation of strcasecmp()
   # do not call AC_REPLACE_FUNCS, as it calls AC_CHECK_FUNCS
   AC_CONFIG_LIBOBJ_DIR([src/libs])
   AC_LIBSOURCE([strcasecmp.c])
   AC_LIBOBJ([strcasecmp])
fi

AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
                    [Enable PnetCDF internal debug mode. This also enables safe mode.
                     @<:@default: disabled@:>@])],
    [debug=${enableval}], [debug=no]
)

dnl malloc memory allocation tracing relies on tdelete and tsearch
AC_CHECK_HEADERS([search.h])
if test "x$ac_cv_header_search_h" = xyes ; then
   AC_CHECK_DECL([tsearch], [], [], [[#include <search.h>]])
   if test "x$ac_cv_have_decl_tsearch" = xyes ; then
      AC_CHECK_FUNCS([tsearch])
   fi
   AC_CHECK_DECL([tdelete], [], [], [[#include <search.h>]])
   if test "x$ac_cv_have_decl_tdelete" = xyes ; then
      AC_CHECK_FUNCS([tdelete])
   fi
fi

PNETCDF_DEBUG=0
if test "x${debug}" = xyes; then
   dnl check required functions for enabling malloc tracing
   if test "x${ac_cv_func_tsearch}" = xyes && test "x${ac_cv_func_tdelete}" = xyes ; then
      AC_DEFINE(PNC_MALLOC_TRACE)
   fi
   PNETCDF_DEBUG=1
fi
AC_SUBST(PNETCDF_DEBUG)
AM_CONDITIONAL(PNETCDF_DEBUG, [test x"$PNETCDF_DEBUG" = x1])

dnl Data type MPI_Offset was first introduced in MPI 2 and since PnetCDF
dnl requires MPI-IO (also introduced in MPI 2), there is no need to check
dnl whether MPI_Offset is defined in mpi.h
dnl AC_CHECK_TYPE([MPI_Offset], [], [], [#include <mpi.h>])

dnl check the size of MPI_Offset. PnetCDF requires it be 8 bytes.
AC_CHECK_SIZEOF([MPI_Offset], [], [#include <mpi.h>])
if test "$ac_cv_sizeof_MPI_Offset" -lt "8"; then
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
     Building of PnetCDF cannot continue due to the size of MPI_Offset
     being less than 8 bytes. Please use an MPI library that supports
     large files. Abort.
   -----------------------------------------------------------------------])
fi

AC_CHECK_SIZEOF([MPI_Aint], [], [#include <mpi.h>])
AM_CONDITIONAL(SIZEOF_MPI_AINT_IS_4, [test x$ac_cv_sizeof_MPI_Aint = x4])

dnl the nonblocking routines build up lists of requests with MPI_Type_struct.
dnl If MPI_Offset not the same size as MPI_Aint, the arrays passed around will
dnl get mangled.
if test "$ac_cv_sizeof_MPI_Offset" -ne "$ac_cv_sizeof_MPI_Aint"; then
   AC_MSG_WARN([
   -----------------------------------------------------------------------
     MPI_Offset and MPI_Aint are detected of different sizes.
     The request aggregation feature implemented in non-blocking APIs is
     thus disabled.
   -----------------------------------------------------------------------])
else
   AC_DEFINE(ENABLE_REQ_AGGREGATION)
   ENABLE_REQ_AGGREGATION=1
   AC_SUBST(ENABLE_REQ_AGGREGATION)
fi

dnl check availability of MPI_Count, which is added in MPI 3.0 standard.
dnl The new MPI functions added in 3.0 that have arguments of type MPI_Count:
dnl MPI_Type_size_x, MPI_Type_get_extent_x, MPI_Type_get_true_extent_x, and
dnl MPI_Get_elements_x
dnl AC_CHECK_TYPE([MPI_Count], [], [], [#include <mpi.h>])
dnl if test "x${ac_cv_type_MPI_Count}" = xyes; then
dnl    AC_DEFINE(HAVE_MPI_COUNT)
dnl fi

dnl Because MPI-IO was first introduced in MPI 2.0 and MPI-IO is a requirement
dnl of PnetCDF, there is no need to check the below APIs which were also
dnl introduced in MPI 2.0.
dnl AC_CHECK_FUNCS([MPI_Get_address \
dnl                 MPI_Type_create_subarray \
dnl                 MPI_Type_create_hvector \
dnl                 MPI_Type_create_hindexed \
dnl                 MPI_Type_create_struct \
dnl                 MPI_Type_create_resized \
dnl                 MPI_Type_get_extent])

dnl Check presence of MPI COMBINERS. These are of type int.
dnl These are introduced in MPI 2.0. As PnetCDF requires an MPI library that
dnl supports MPI-IO and MPI-IO was first introduced in MPI 2.0, checking these
dnl MPI constants is not necessary.
dnl UD_CHECK_MPI_CONSTANTS([MPI_COMBINER_NAMED,
dnl                         MPI_COMBINER_DUP,
dnl                         MPI_COMBINER_CONTIGUOUS,
dnl                         MPI_COMBINER_VECTOR,
dnl                         MPI_COMBINER_HVECTOR,
dnl                         MPI_COMBINER_HINDEXED,
dnl                         MPI_COMBINER_INDEXED_BLOCK,
dnl                         MPI_COMBINER_STRUCT,
dnl                         MPI_COMBINER_SUBARRAY,
dnl                         MPI_COMBINER_DARRAY,
dnl                         MPI_COMBINER_F90_REAL,
dnl                         MPI_COMBINER_F90_COMPLEX,
dnl                         MPI_COMBINER_F90_INTEGER,
dnl                         MPI_COMBINER_RESIZED,
dnl                         MPI_COMBINER_HINDEXED_BLOCK],
dnl                         [], [], [[#include <mpi.h>]])

dnl Below MPI constants have been deprecated since MPI 3.0
dnl Note from autoconfig user guide: "Unlike the other AC_CHECK_*S macros, when
dnl a symbol is not declared, HAVE_DECL_symbol is defined to 0 instead of
dnl leaving HAVE_DECL_symbol undeclared." So, use it like this example:
dnl #if defined MPI_COMBINER_NAMED && MPI_COMBINER_NAMED
dnl     ... use MPI_COMBINER_NAMED
dnl #endif
dnl However, we cannot use AC_CHECK_DECLS to check whether these MPI constants
dnl are deprecated. OpenMPI still define those constants deprecated by MPI 3.0
dnl as error messages. We need to call AC_COMPILE_IFELSE to check whether they
dnl can be used without compilation errors.
dnl
dnl Check MPI constants below only when MPI library implements MPI 3 and later
dnl as they were introduced in MPI 2, but deprecated since MPI 3.
if test "$mpi_version" -ge "3" ; then
   UD_CHECK_MPI_CONSTANTS([MPI_COMBINER_HVECTOR_INTEGER,
                           MPI_COMBINER_HINDEXED_INTEGER,
                           MPI_COMBINER_STRUCT_INTEGER],
                          [], [], [[#include <mpi.h>]])
fi

dnl Check presence of various MPI error classes. Introduced in MPI 2.0.
dnl These could be enums, so we have to do compile checks.
dnl AC_CHECK_DECLS([MPI_ERR_FILE_EXISTS,
dnl                 MPI_ERR_NO_SUCH_FILE,
dnl                 MPI_ERR_AMODE,
dnl                 MPI_ERR_NOT_SAME,
dnl                 MPI_ERR_BAD_FILE,
dnl                 MPI_ERR_READ_ONLY,
dnl                 MPI_ERR_ACCESS,
dnl                 MPI_ERR_NO_SPACE,
dnl                 MPI_ERR_QUOTA],
dnl                 [], [], [[#include <mpi.h>]])

dnl Check presence of C MPI data types.
dnl Note they are of type MPI_Datatype, thus cannot use AC_CHECK_DECLS
dnl UD_CHECK_MPI_DATATYPE(MPI_CHAR)		dnl first defined in MPI 1.0
dnl UD_CHECK_MPI_DATATYPE(MPI_BYTE)		dnl first defined in MPI 1.0
dnl UD_CHECK_MPI_DATATYPE(MPI_SIGNED_CHAR)	dnl first defined in MPI 2.0
dnl UD_CHECK_MPI_DATATYPE(MPI_UNSIGNED_CHAR)	dnl first defined in MPI 1.0
dnl UD_CHECK_MPI_DATATYPE(MPI_SHORT)		dnl first defined in MPI 1.0
dnl UD_CHECK_MPI_DATATYPE(MPI_UNSIGNED_SHORT)	dnl first defined in MPI 1.0
dnl UD_CHECK_MPI_DATATYPE(MPI_INT)		dnl first defined in MPI 1.0
dnl UD_CHECK_MPI_DATATYPE(MPI_UNSIGNED)		dnl first defined in MPI 1.0
dnl UD_CHECK_MPI_DATATYPE(MPI_LONG)		dnl first defined in MPI 1.0
dnl UD_CHECK_MPI_DATATYPE(MPI_UNSIGNED_LONG)	dnl first defined in MPI 1.0
dnl UD_CHECK_MPI_DATATYPE(MPI_FLOAT)		dnl first defined in MPI 1.0
dnl UD_CHECK_MPI_DATATYPE(MPI_DOUBLE)		dnl first defined in MPI 1.0
dnl UD_CHECK_MPI_DATATYPE(MPI_LONG_LONG_INT)	dnl first defined in MPI 1.0
dnl UD_CHECK_MPI_DATATYPE(MPI_UNSIGNED_LONG_LONG) dnl first defined in MPI 2.0
dnl UD_CHECK_MPI_DATATYPE(MPI_UB)		dnl first defined in MPI 1.0
dnl UD_CHECK_MPI_DATATYPE(MPI_LB)		dnl first defined in MPI 1.0

dnl C data type MPI_Offset was first introduced in MPI 2.0, but MPI_OFFSET was
dnl introduced as a predefined MPI datatype in MPI 2.2
if test "$mpi_version" -ge "3" ; then
   ac_cv_CHECK_MPI_DATATYPE_MPI_OFFSET=yes
else
   UD_CHECK_MPI_DATATYPE(MPI_OFFSET)
fi
AM_CONDITIONAL(DECL_MPI_OFFSET, [test x$ac_cv_CHECK_MPI_DATATYPE_MPI_OFFSET = xyes])

dnl Check for presence of Fortran types
dnl These could be enums, so we have to do compile checks.
dnl
dnl We do this for a couple of reasons.  First, the MPI might have been
dnl built without Fortran support, in which case these types might not
dnl exist.  Second, we need to map these types to corresponding C types
dnl where possible to simplify processing at run time.
dnl if test "x${has_fortran}" = xyes ; then
dnl     UD_CHECK_MPI_DATATYPE(MPI_CHARACTER)		dnl first defined in MPI 1.0
dnl     UD_CHECK_MPI_DATATYPE(MPI_INTEGER)		dnl first defined in MPI 1.0
dnl     UD_CHECK_MPI_DATATYPE(MPI_INTEGER1)		dnl first defined in MPI 1.0
dnl     UD_CHECK_MPI_DATATYPE(MPI_INTEGER2)		dnl first defined in MPI 1.0
dnl     UD_CHECK_MPI_DATATYPE(MPI_INTEGER4)		dnl first defined in MPI 1.0
dnl     UD_CHECK_MPI_DATATYPE(MPI_INTEGER8)		dnl first defined in MPI 2.1
dnl     UD_CHECK_MPI_DATATYPE(MPI_REAL)			dnl first defined in MPI 1.0
dnl     UD_CHECK_MPI_DATATYPE(MPI_REAL2)		dnl first defined in MPI 1.0
dnl     UD_CHECK_MPI_DATATYPE(MPI_REAL4)		dnl first defined in MPI 1.0
dnl     UD_CHECK_MPI_DATATYPE(MPI_REAL8)		dnl first defined in MPI 1.0
dnl     UD_CHECK_MPI_DATATYPE(MPI_DOUBLE_PRECISION)	dnl first defined in MPI 1.0
dnl fi

AC_C_CHAR_UNSIGNED
AC_C_BIGENDIAN
AM_CONDITIONAL(IS_BIGENDIAN, [test x$ac_cv_c_bigendian = xyes])
AC_SUBST(ac_cv_c_bigendian)dnl for src/utils/pnetcdf-config.in

AC_ARG_ENABLE([in-place-swap],
    [AS_HELP_STRING([--(en/dis)able-in-place-swap],
                    [Enable/disable memory in-place byte swap on Little Endian
                     machines. @<:@default: auto@:>@])],
    [in_place_swap=${enableval}], [in_place_swap=auto]
)
UD_MSG_DEBUG([in_place_swap=$in_place_swap])
dnl ENABLE_IN_PLACE_SWAP=0
dnl if test "x${in_place_swap}" = xyes ; then
dnl    ENABLE_IN_PLACE_SWAP=1
dnl fi
dnl AC_DEFINE(ENABLE_IN_PLACE_SWAP, $ENABLE_IN_PLACE_SWAP)
dnl DISABLE_IN_PLACE_SWAP=0
dnl if test "x${in_place_swap}" = xno ; then
dnl    if test "x${ac_cv_c_bigendian}" = xyes ; then
dnl       AC_MSG_WARN([--disable-in-place-swap takes no effect on Big Endian])
dnl       in_place_swap=yes
dnl    else
dnl       DISABLE_IN_PLACE_SWAP=1
dnl    fi
dnl fi
dnl AC_DEFINE(DISABLE_IN_PLACE_SWAP, $DISABLE_IN_PLACE_SWAP)

dnl IN_PLACE_SWAP being -1 corresponds to "auto"
IN_PLACE_SWAP=-1
if test "x${ac_cv_c_bigendian}" = xyes && test "x${in_place_swap}" != xauto ; then
   AC_MSG_WARN([--(en/dis)able-in-place-swap takes no effect on Big Endian])
else
   if test "x${in_place_swap}" = xyes ; then
      dnl IN_PLACE_SWAP being 1 corresponds to "yes"
      IN_PLACE_SWAP=1
   fi
   if test "x${in_place_swap}" = xno ; then
      dnl IN_PLACE_SWAP being 0 corresponds to "no"
      IN_PLACE_SWAP=0
   fi
fi
AC_SUBST(IN_PLACE_SWAP, $IN_PLACE_SWAP)dnl for src/include/pnetcdf-h.in
dnl in_place_swap can be yes, no, auto
AC_SUBST(in_place_swap)dnl for src/utils/pnetcdf-config.in

dnl For big Endian, put buffer needs no byte swap and hence can be declared as
dnl INTENT(IN). For little Endian, put buffer may be used for byte swap in
dnl place and hence must be declared as INTENT(INOUT).
dnl This will configure/produce the file src/binding/f90/api.f90
if test "x${ac_cv_c_bigendian}" = xyes || test "x${in_place_swap}" = xno ; then
   INTENTV="IN"
else
   INTENTV="INOUT"
fi
AC_SUBST(INTENTV)dnl for src/binding/f90/api.fh.in

AC_TYPE_SIZE_T
AC_TYPE_OFF_T
AC_TYPE_SSIZE_T

AC_CHECK_TYPES([ptrdiff_t, schar, uchar, ushort, uint, longlong, ulonglong, int64, uint64])

AC_CHECK_SIZEOF(size_t)
AC_CHECK_SIZEOF(off_t)
AC_CHECK_SIZEOF(char)
AC_CHECK_SIZEOF(signed char)
AC_CHECK_SIZEOF(unsigned char)
AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(unsigned short int)
AC_CHECK_SIZEOF(unsigned short)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(long long)
AC_CHECK_SIZEOF(unsigned long long)
AC_CHECK_SIZEOF(ptrdiff_t)

if test "$ac_cv_type_schar" = yes ; then
   AC_CHECK_SIZEOF(schar)
fi
if test "$ac_cv_type_uchar" = yes ; then
   AC_CHECK_SIZEOF(uchar)
fi
if test "$ac_cv_type_ushort" = yes ; then
   AC_CHECK_SIZEOF(ushort)
fi
if test "$ac_cv_type_uint" = yes ; then
   AC_CHECK_SIZEOF(uint)
fi
if test "$ac_cv_type_longlong" = yes ; then
   AC_CHECK_SIZEOF(longlong)
fi
if test "$ac_cv_type_ulonglong" = yes ; then
   AC_CHECK_SIZEOF(ulonglong)
fi

if test "x${has_fortran}" = xyes ; then
   dnl check Fortran default integer size
   dnl PnetCDF currently does not support the default integer of size 8
   AC_CHECK_SIZEOF([MPI_Fint], [], [#include <mpi.h>])
   UD_MSG_DEBUG(["sizeof(MPI_Fint) = $ac_cv_sizeof_MPI_Fint"])
   if test "$ac_cv_sizeof_MPI_Fint" = 8 ; then
      AC_MSG_ERROR([PnetCDF does not support Fortran default integer size of 8 bytes])
   fi

   if test "$cross_compiling" = yes; then
      UD_CHECK_FORTRAN_TYPE([NF_INT1_T], [integer*1 byte "integer(kind=1)"])
      UD_CHECK_FORTRAN_TYPE([NF_INT2_T], [integer*2 "integer(kind=2)"])
      UD_CHECK_FORTRAN_TYPE([NF_INT8_T], [integer*8 "integer(kind=8)"])
   else
      UD_FORTRAN_TYPES
   fi
   AC_SUBST(NF_INT1_T)
   AC_SUBST(NF_INT2_T)
   AC_SUBST(NF_INT8_T)

   dnl NFMPI_OFFSET="integer*$ac_cv_sizeof_MPI_Offset"
   dnl AC_MSG_CHECKING([for Fortran NFMPI_OFFSET "$NFMPI_OFFSET"])
   dnl dnl "
   dnl AC_LANG_PUSH([Fortran 77])
   dnl AC_COMPILE_IFELSE(
   dnl    [AC_LANG_SOURCE([
   dnl          subroutine sub(value)
   dnl          $NFMPI_OFFSET value
   dnl          end
   dnl    ])],
   dnl    [ac_cv_NFMPI_OFFSET=yes], [ac_cv_NFMPI_OFFSET=no]
   dnl )
   dnl AC_LANG_POP([Fortran 77])
   dnl if test "$ac_cv_NFMPI_OFFSET" = yes ; then
   dnl    AC_MSG_RESULT(yes)
   dnl else
   dnl    AC_MSG_RESULT(no)
   dnl    AC_MSG_ERROR([F77 does not support "$NFMPI_OFFSET"])
   dnl    dnl "
   dnl fi
   dnl ${RM} -rf conftest*
fi

SIZEOF_MPI_OFFSET=$ac_cv_sizeof_MPI_Offset
AC_SUBST(SIZEOF_MPI_OFFSET)

HAVE_F77_GNU_INT=no
HAVE_F77_INT1=no
HAVE_F77_INT2=no
HAVE_F77_INT8=no
if test "x${has_fortran}" = xyes ; then
   UD_CHECK_F77_GNU_INT
   UD_MSG_DEBUG([ac_cv_f77_gnu_int=$ac_cv_f77_gnu_int])
   HAVE_F77_GNU_INT=$ac_cv_f77_gnu_int
   AC_SUBST(HAVE_F77_GNU_INT)

   UD_CHECK_F77_INT1
   UD_MSG_DEBUG([ac_cv_f77_int1=$ac_cv_f77_int1])
   HAVE_F77_INT1=$ac_cv_f77_int1
   AC_SUBST(HAVE_F77_INT1)
   UD_CHECK_F77_INT2
   UD_MSG_DEBUG([ac_cv_f77_int2=$ac_cv_f77_int2])
   HAVE_F77_INT2=$ac_cv_f77_int2
   AC_SUBST(HAVE_F77_INT2)
   UD_CHECK_F77_INT8
   UD_MSG_DEBUG([ac_cv_f77_int8=$ac_cv_f77_int8])
   HAVE_F77_INT8=$ac_cv_f77_int8
   AC_SUBST(HAVE_F77_INT8)
fi
AM_CONDITIONAL(HAVE_F77_GNU_INT, [test x$HAVE_F77_GNU_INT = xyes])
AM_CONDITIONAL(HAVE_F77_INT1, [test x$HAVE_F77_INT1 = xyes])
AM_CONDITIONAL(HAVE_F77_INT2, [test x$HAVE_F77_INT2 = xyes])
AM_CONDITIONAL(HAVE_F77_INT8, [test x$HAVE_F77_INT8 = xyes])

dnl UD_MAKEWHATIS

dnl
dnl GNU coverage
dnl
dnl This is for internal testing only. It should not be enabled for building a
dnl production PnetCDF. This is because running an executable compiled with
dnl coverage will produce an output file named "gmon.out". Since coverage is
dnl not parallelized, running a program compiled with coverage may cause
dnl problems on concurrently writing to gmon.out in conflicts, possible
dnl corrupting the file or program hanging. Thus, make target "ptest" should
dnl also be disabled when coverage is enabled.
dnl
dnl After all other tests, optionally enable coverage,  we do this last
dnl because legend has it that sometimes on some compilers the coverage flags
dnl mess up other checks
dnl
if test "x${has_fortran}" = xyes ; then
   enable_f77=yes
   enable_fc=yes
fi
PAC_ENABLE_COVERAGE

dnl Some C++ compilers (e.g. Fujitsu) have macro __FUNCTION__, but not __func__
use_cxx_macro_function=no
if test "x${has_mpicxx}" = xyes ; then
   UD_CXX_MACRO_FUNC
   if test "x${ac_cv_cxx_macro_func}" = xno && test "x${ac_cv_cxx_macro_function}" = xyes ; then
      use_cxx_macro_function=yes
   fi
fi
AM_CONDITIONAL(REPLACE_CXX_FUNC_MACRO, [test x$use_cxx_macro_function = xyes])

if test "x${debug}" = xyes; then
   dnl add -g flag if not presented
   dnl remove all -O and -fast flags
   dnl add -O0 to all flags
   # check exit status of grep command is more portable than using -q
   str_found=`echo "${CFLAGS}" | ${EGREP} -- "-g"`
   if test "x$?" != x0 ; then
      CFLAGS="$CFLAGS -g"
   fi
   CFLAGS=`echo $CFLAGS | ${SED} 's/-O. *//g' | ${SED} 's/-fast *//g'`
   CFLAGS="$CFLAGS -O0"

   if test "x${has_mpicxx}" = xyes ; then
      str_found=`echo "${CXXFLAGS}" | ${EGREP} -- "-g"`
      if test "x$?" != x0 ; then
         CXXFLAGS="$CXXFLAGS -g"
      fi
      CXXFLAGS=`echo $CXXFLAGS | ${SED} 's/-O. *//g' | ${SED} 's/-fast *//g'`
      CXXFLAGS="$CXXFLAGS -O0"
   fi

   if test "x${has_fortran}" = xyes ; then
      str_found=`echo "${FFLAGS}" | ${EGREP} -- "-g"`
      if test "x$?" != x0 ; then
         FFLAGS="$FFLAGS -g"
      fi
      str_found=`echo "${FCFLAGS}" | ${EGREP} -- "-g"`
      if test "x$?" != x0 ; then
         FCFLAGS="$FCFLAGS -g"
      fi

       FFLAGS=`echo $FFLAGS  | ${SED} 's/-O. *//g' | ${SED} 's/-fast *//g'`
      FCFLAGS=`echo $FCFLAGS | ${SED} 's/-O. *//g' | ${SED} 's/-fast *//g'`
       FFLAGS="$FFLAGS -O0"
      FCFLAGS="$FCFLAGS -O0"
   fi
   unset str_found
fi

dnl Starting from GNU Fortran 10.0.0, function/subroutine argument type
dnl mismatch becomes a compile error. "Mismatches between actual and dummy
dnl argument lists in a single file are now rejected with an error. Use the new
dnl option -fallow-argument-mismatch to turn these errors into warnings; this
dnl option is implied with -std=legacy. -Wargument-mismatch has been removed."
dnl See https://gcc.gnu.org/gcc-10/changes.html and
dnl https://github.com/Parallel-NetCDF/PnetCDF/issues/61
if test "x${has_fortran}" = xyes && test "x${ax_cv_fc_compiler_vendor}" = xgnu ; then
   AC_MSG_CHECKING([whether FFLAGS already includes -fallow-argument-mismatch])
   str_found=`echo "${FFLAGS}" | ${EGREP} -- "-fallow-argument-mismatch"`
   if test "x$?" = x0 ; then
      has_fallow_argument_mismatch=yes
   else
      has_fallow_argument_mismatch=no
   fi
   AC_MSG_RESULT([$has_fallow_argument_mismatch])

   if test "x${has_fallow_argument_mismatch}" = xno ; then
      AC_MSG_CHECKING([whether MPIF77 accepts option -fallow-argument-mismatch])
      AC_LANG_PUSH([Fortran 77])
      saved_FFLAGS=$FFLAGS
      FFLAGS="$FFLAGS -fallow-argument-mismatch"
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [])],
                        [fallow_argument_mismatch=yes],
                        [fallow_argument_mismatch=no])
      AC_MSG_RESULT([$fallow_argument_mismatch])
      if test "x${fallow_argument_mismatch}" = xno ; then
         FFLAGS=$saved_FFLAGS
      fi
      AC_LANG_POP([Fortran 77])
   fi

   AC_MSG_CHECKING([whether FCFLAGS already includes -fallow-argument-mismatch])
   str_found=`echo "${FCFLAGS}" | ${EGREP} -- "-fallow-argument-mismatch"`
   if test "x$?" = x0 ; then
      has_fallow_argument_mismatch=yes
   else
      has_fallow_argument_mismatch=no
   fi
   AC_MSG_RESULT([$has_fallow_argument_mismatch])

   if test "x${has_fallow_argument_mismatch}" = xno ; then
      AC_MSG_CHECKING([whether MPIF90 accepts option -fallow-argument-mismatch])
      AC_LANG_PUSH([Fortran])
      saved_FCFLAGS=$FCFLAGS
      FCFLAGS="$FCFLAGS -fallow-argument-mismatch"
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [])],
                        [fallow_argument_mismatch=yes],
                        [fallow_argument_mismatch=no])
      AC_MSG_RESULT([$fallow_argument_mismatch])
      if test "x${fallow_argument_mismatch}" = xno ; then
         FCFLAGS=$saved_FCFLAGS
      fi
      AC_LANG_POP([Fortran])
   fi
fi

chmod u+x ${srcdir}/scripts/install-sh

AC_ARG_ENABLE([subfiling],
   [AS_HELP_STRING([--enable-subfiling],
                   [Enable subfiling support. @<:@default: disabled@:>@])],
   [enable_subfiling=${enableval}], [enable_subfiling=no]
)

ENABLE_SUBFILING=0
if test "x$enable_subfiling" = "xyes" ; then
   AC_DEFINE(ENABLE_SUBFILING)
   ENABLE_SUBFILING=1
fi
AC_SUBST(ENABLE_SUBFILING)
AM_CONDITIONAL(ENABLE_SUBFILING, [test x$enable_subfiling = xyes])

AC_ARG_ENABLE([thread-safe],
    [AS_HELP_STRING([--enable-thread-safe],
                    [Enable thread-safe capability. @<:@default: disabled@:>@])],
    [thread_safe=${enableval}], # --enable-thread-safe or --disable-thread-safe is used
    [thread_safe=no]            # neither --enable-thread-safe nor --disable-thread-safe is used
)
ENABLE_THREAD_SAFE=0
if test "x${thread_safe}" = xyes ; then
   AC_DEFINE(ENABLE_THREAD_SAFE)
   ENABLE_THREAD_SAFE=1
fi
AC_SUBST(ENABLE_THREAD_SAFE)
AM_CONDITIONAL(ENABLE_THREAD_SAFE, [test x$thread_safe = xyes])

AC_ARG_WITH([pthread],
   [AS_HELP_STRING([--with-pthread=DIR],
                   [Search Pthreads library within the supplied path DIR,
                    when --enable-thread-safe is enabled.])],
   [ dnl this clause is run when --with-pthread or --without-pthread is used
   if test "x$withval" = xno && test "x$thread_safe" = xyes ; then
      AC_MSG_ERROR([conflicted options used: --enable-thread-safe and --without-pthread])
   elif test "x$withval" = xyes ; then
      pthread_dir=check
   else
      pthread_dir=$withval
   fi
   ],[ dnl this clause is run when neither --with-pthread nor --without-pthread is used
   pthread_dir=check
   ]
)

if test "x${thread_safe}" = xyes ; then
   if test "x${pthread_dir}" = xcheck ; then
      AC_CHECK_HEADER([pthread.h],,[AC_MSG_ERROR([Cannot find pthread.h])])
      AC_CHECK_LIB([pthread], [pthread_self],, [AC_MSG_ERROR([Cannot find pthread library])])
   else
      case "$pthread_dir" in
         *,*)
            pthread_inc="`echo $pthread_dir | cut -f1 -d,`"
            pthread_lib="`echo $pthread_dir | cut -f2 -d, -s`"
            ;;
         *)
            # check $pthread_dir for empty is done above
            pthread_inc="$pthread_dir/include"
            pthread_lib="$pthread_dir/lib"
            ;;
      esac
      CPPFLAGS="$CPPFLAGS -I$pthread_inc"
      AC_CHECK_HEADER([pthread.h],, [AC_MSG_ERROR([Cannot find pthread.h under $pthread_dir])])
      LDFLAGS="$LDFLAGS -L$pthread_lib"
      AC_CHECK_LIB([pthread], [pthread_self],, [AC_MSG_ERROR([Cannot find pthread library under $pthread_dir])])
   fi
fi

enable_netcdf4=no
NETCDF4_INSTALL=""
netcdf_includedir=""
netcdf4_libdir=""
AC_ARG_WITH(netcdf4,
   [AS_HELP_STRING([--with-netcdf4@<:@=INC,LIB | =DIR@:>@],
                   [Enable NetCDF-4 feature and provide the NetCDF-4 installation path(s):
    --with-netcdf4=INC,LIB for include and lib paths separated by a comma.
    --with-netcdf4=DIR for the path containing include/ and lib/ subdirectories.
    @<:@default: disabled@:>@
   ])],[ dnl this clause is run when --with-netcdf4 or --without-netcdf4 is used
   if test "x${withval}" != xyes && test "x${withval}" != xno ; then
      NETCDF4_INSTALL=${withval}
      enable_netcdf4=yes

      case $NETCDF4_INSTALL in
        "") ;;
        *,*)
           netcdf_includedir="`echo $NETCDF4_INSTALL |cut -f1 -d,`"
           netcdf4_libdir="`echo $NETCDF4_INSTALL |cut -f2 -d, -s`"
           netcdf_bindir="$netcdf_includedir/../bin"
           ;;
        *)
           netcdf_includedir="$NETCDF4_INSTALL/include"
           netcdf4_libdir="$NETCDF4_INSTALL/lib"
           netcdf_bindir="$NETCDF4_INSTALL/bin"
           ;;
      esac
   elif test "x${withval}" = xyes ; then
      dnl user may add NetCDF INC and LIB to CFLAGS and LDFLAGS, respectively.
      enable_netcdf4=yes
   fi
   UD_MSG_DEBUG([NETCDF4_INSTALL=$NETCDF4_INSTALL])
   UD_MSG_DEBUG([netcdf_includedir=$netcdf_includedir])
   UD_MSG_DEBUG([netcdf4_libdir=$netcdf4_libdir])
   ]
)

if test "x$enable_netcdf4" = xyes ; then
   # check availability of nc-config only under $NETCDF4_INSTALL/bin
   nc_config=
   if test "x$NETCDF4_INSTALL" != x ; then
      AC_PATH_PROG([nc_config],[nc-config],,[$netcdf_bindir])
      if test "x$nc_config" = x ; then
         # nc-config may not be installed, even if NetCDF4 is installed. For
         # example, there are only include and lib installed under
         # NETCDF4_INSTALL on Cori @NERSC. Missing bin folder under
         # NETCDF4_INSTALL may be due to the cross compile environment to make
         # all executables, such as ncdump and ncgen, supposed to be under bin
         # only available for login nodes and thus those executables may be
         # installed in a different location. In this case, we just use and
         # check header and library files under $NETCDF4_INSTALL.

         UD_MSG_DEBUG([nc-config is not available under $NETCDF4_INSTALL/bin])
         dnl AC_MSG_WARN([
         dnl ------------------------------------------------------------
         dnl Cannot find nc_config, a utility command to show the
         dnl configuration of the NetCDF-4 library. It is required to
         dnl build PnetCDF with NetCDF-4 support. Use option
         dnl     --with-netcdf4=/path/to/implementation
         dnl to specify the location of NetCDF-4 build. In addition,
         dnl please make sure the MPI C compiler is compatible with the
         dnl one used to build NetCDF-4. The NetCDF-4 library used must
         dnl be built with parallel I/O and NETCDF-4 support. Check
         dnl 'config.log' for more information. Stopping ...
         dnl ------------------------------------------------------------])
      fi
   else
      dnl Check nc-config under $PATH
      AC_PATH_PROG([nc_config],[nc-config])
      if test "x$nc_config" != x ; then
         netcdf_includedir=`$nc_config --includedir`
         UD_MSG_DEBUG([netcdf_includedir=$netcdf_includedir])
         netcdf4_libdir=`$nc_config --libdir`
         UD_MSG_DEBUG([netcdf4_libdir=$netcdf4_libdir])
      fi
   fi

   if test "x$nc_config" != x ; then
      dnl Check if NetCDF is built with HDF5 enabled
      netcdf_has_hdf5=`$nc_config --has-hdf5`
      if test "x$netcdf_has_hdf5" = xno ; then
         AC_MSG_ERROR([
         ------------------------------------------------------------
         NetCDF-4 library is not built with HDF5 enabled, as indicated
         by command '$nc_config --has-hdf5'. Stop
         ------------------------------------------------------------])
      fi
      unset netcdf_has_hdf5
   fi

   dnl find out how NetCDF is configured.

   dnl No need to add the include paths used to build NetCDF to CPPFLAGS
   dnl permanently, as it is not needed to compile PnetCDF codes
   saved_CPPFLAGS=$CPPFLAGS
   if test "x$netcdf_includedir" != x ; then
       CPPFLAGS="-I${netcdf_includedir} ${CPPFLAGS}"
   fi
   have_netcdf4=yes
   AC_CHECK_HEADER([netcdf.h],      [], [have_netcdf4=no])
   AC_CHECK_HEADERS([netcdf_meta.h],[], [have_netcdf4=no], [[#include <netcdf.h>]])
   if test "x${have_netcdf4}" = xno; then
      AC_MSG_ERROR([
      ------------------------------------------------------------
      Missing NetCDF-4 header files 'netcdf.h' or 'netcdf_meta.h'
      required to build PnetCDF with NetCDF-4 support. Use option
          --with-netcdf4=/path/to/implementation
      to specify the location of NetCDF-4 installation. In addition,
      please make sure the MPI C compiler is compatible with the one
      used to build NetCDF-4. Check 'config.log' for more information.
      Stopping ...
      ------------------------------------------------------------])
   fi
   AC_CHECK_HEADERS([netcdf_par.h],  [], [have_netcdf4=no], [[#include <netcdf.h>]])
   if test "x${have_netcdf4}" = xno; then
      AC_MSG_ERROR([
      ------------------------------------------------------------
      Missing NetCDF-4 header file 'netcdf_par.h' required by PnetCDF.
      The NetCDF-4 library must be built with parallel I/O and
      NETCDF-4 support. Check 'config.log' for more information.
      Stopping ...
      ------------------------------------------------------------])
   fi

   dnl Check if NetCDF version is 4.6.2 or later
   AC_MSG_CHECKING([whether NetCDF version is 4.6.2 or later])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <netcdf_meta.h>
#define VER_NUM (NC_VERSION_MAJOR*1000000 + NC_VERSION_MINOR*1000 + NC_VERSION_PATCH)
#if (VER_NUM < 4006002)
#error NetCDF version is older than 4.6.2
#endif
   ]])], [netcdf_ge_4_6_2=yes], [netcdf_ge_4_6_2=no])
   AC_MSG_RESULT([$netcdf_ge_4_6_2])

   dnl Check if NetCDF version is 4.8.0 or 4.8.1
   AC_MSG_CHECKING([whether NetCDF version is 4.8.0 or 4.8.1])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <netcdf_meta.h>
#define VER_NUM (NC_VERSION_MAJOR*1000000 + NC_VERSION_MINOR*1000 + NC_VERSION_PATCH)
#if (VER_NUM != 4008000) && (VER_NUM != 4008001)
#error NetCDF version is neither 4.8.0 nor 4.8.1
#endif
   ]])], [netcdf_480_481=yes], [netcdf_480_481=no])
   AC_MSG_RESULT([$netcdf_480_481])

   dnl Check if NetCDF is built with parallel I/O enabled
   AC_MSG_CHECKING([whether NetCDF library is built with parallel I/O enabled])
   if test x$netcdf_ge_4_6_2 = xyes; then
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <netcdf_meta.h>
#if !defined(NC_HAS_PARALLEL4) || NC_HAS_PARALLEL4 == 0
#error NetCDF NC_HAS_PARALLEL4 is either not defined or defined to 0
#endif
      ]])], [netcdf_has_parallel4=yes], [netcdf_has_parallel4=no])
      if test "x$netcdf_has_parallel4" = xno ; then
         AC_MSG_RESULT([no])
         AC_MSG_ERROR([NetCDF-4 library is not built with HDF5 parallel I/O enabled])
      fi
   else
      dnl NetCDF 4.6.1 and earlier do not have --has-parallel4 configure option.
      dnl The test below checking whether HDF5 parallel feature is enabled will
      dnl miss the case when NetCDF is built with parallel HDF5 disabled and
      dnl PnetCDF enabled, i.e. configured with options --enable-pnetcdf
      dnl --enable-netcdf-4 --disable-parallel4, which will make command
      dnl "nc-config --has-parallel" return yes.
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <netcdf_meta.h>
#if !defined(NC_HAS_PARALLEL) || NC_HAS_PARALLEL == 0
#error NetCDF NC_HAS_PARALLEL is either not defined or defined to 0
#endif
      ]])], [netcdf_has_parallel=yes], [netcdf_has_parallel=no])
      if test "x$netcdf_has_parallel" = xno ; then
         AC_MSG_RESULT([no])
         AC_MSG_ERROR([NetCDF-4 library is not built with parallel I/O enabled])
      fi
   fi
   AC_MSG_RESULT([yes])

   dnl Check if NetCDF version is 4.5.0 or later
   AC_MSG_CHECKING([whether NetCDF version is 4.5.0 or later])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <netcdf_meta.h>
#if (NC_VERSION_MAJOR*1000000 + NC_VERSION_MINOR*1000 + NC_VERSION_PATCH < 4005000)
#error NetCDF version is older than 4.5.0
#endif
   ]])], [netcdf_ge_4_5_0=yes], [netcdf_ge_4_5_0=no])
   AC_MSG_RESULT([$netcdf_ge_4_5_0])
   if test "x$netcdf_ge_4_5_0" = xyes ; then
       AC_DEFINE(NETCDF_GE_4_5_0)
   fi

   dnl Check if NetCDF version is 4.4.2 or later
   dnl nc-config option --libdir is first added in NetCDF 4.4.2
   AC_MSG_CHECKING([whether NetCDF version is 4.4.2 or later])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <netcdf_meta.h>
#if (NC_VERSION_MAJOR*1000000 + NC_VERSION_MINOR*1000 + NC_VERSION_PATCH < 4004002)
#error NetCDF version is older than 4.4.2
#endif
   ]])], [netcdf_ge_4_4_2=yes], [netcdf_ge_4_4_2=no])
   AC_MSG_RESULT([$netcdf_ge_4_4_2])

   AC_MSG_CHECKING([whether relax-coord-bound is enabled in NetCDF-4])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <netcdf_meta.h>
#if !defined(NC_RELAX_COORD_BOUND) || NC_RELAX_COORD_BOUND == 0
#error NetCDF is built with relax-coord-bound disable
#endif
   ]])], [nc_relax_coord_bound=yes], [nc_relax_coord_bound=no])
   AC_MSG_RESULT([$nc_relax_coord_bound])

   # Build on top of NetCDF4 with PnetCDF already enabled is not supported
   AC_MSG_CHECKING([whether NetCDF is built with PnetCDF])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <netcdf_meta.h>
#if !defined(NC_HAS_PNETCDF) || NC_HAS_PNETCDF == 0
#error NetCDF is not built with PnetCDF enabled
#endif
   ]])], [netcdf_has_pnetcdf=yes], [netcdf_has_pnetcdf=no])
   AC_MSG_RESULT([$netcdf_has_pnetcdf])
   if test x$netcdf_has_pnetcdf = xyes; then
      AC_MSG_ERROR([NetCDF-4 built with PnetCDF enabled is not supported])
   fi

   if test "x$nc_config" != x ; then
      # Get installation directory path of netcdf library files
      # if test "x$netcdf_ge_4_4_2" = xyes ; then
      #    # Option --libdir is first added in NetCDF 4.4.2
      #    netcdf4_libdir=`$nc_config --libdir`
      #    UD_MSG_DEBUG([netcdf4_libdir=$netcdf4_libdir])
      # fi

      # Get extra libraries used by netcdf library (dependent libraries)
      netcdf_libs=`$nc_config --libs`
      UD_MSG_DEBUG([netcdf_libs=$netcdf_libs])
      dnl remove -lpnetcdf if it contains it
      netcdf_libs=`echo $netcdf_libs | $SED "s/-lpnetcdf//g"`
   else
      # if test "x$NETCDF4_INSTALL" != x ; then
      #    netcdf4_libdir="$NETCDF4_INSTALL/lib"
      #    UD_MSG_DEBUG([netcdf4_libdir=$netcdf4_libdir])
      # fi
      netcdf_libs="-lnetcdf"
      UD_MSG_DEBUG([netcdf_libs=$netcdf_libs])
   fi

   # Validate netcdf library
   saved_LDFLAGS=$LDFLAGS
   saved_LIBS=$LIBS
   if test "x$netcdf4_libdir" != x ; then
      LDFLAGS="-L$netcdf4_libdir $LDFLAGS"
   fi
   LIBS="$netcdf_libs $LIBS"
   AC_SEARCH_LIBS([nc_open], [netcdf],,[AC_MSG_ERROR([Cannot find nc_open in NetCDF-4 library])])
   LDFLAGS=$saved_LDFLAGS
   LIBS=$saved_LIBS
   # Do not add -L$netcdf4_libdir to LDFLAGS and -lnetcdf to LIBS
   # They are only needed when linking executables.

   # Get directory path of netcdf library files
   # netcdf_lib=`$nc_config --libs`
   # UD_MSG_DEBUG([netcdf_lib=$netcdf_lib])
   # dnl separate -L and -l
   # netcdf_ldflags=`echo $netcdf_lib | $SED "s/-l[[^ ]]*//g"`
   # netcdf_libs=`echo $netcdf_lib | $SED "s/-L[[^ ]]*//g"`
   # dnl remove -lpnetcdf if it contains it
   # netcdf_libs=`echo $netcdf_libs | $SED "s/-lpnetcdf//g"`
   # dnl remove extra whitespace
   # netcdf_ldflags=`echo $netcdf_ldflags | $SED "s/[[ ]]* / /g"`
   # netcdf_libs=`echo $netcdf_libs | $SED "s/[[ ]]* / /g"`
   # UD_MSG_DEBUG([netcdf_ldflags=$netcdf_ldflags])
   # UD_MSG_DEBUG([netcdf_libs=$netcdf_libs])
   # LDFLAGS="$LDFLAGS ${netcdf_ldflags}"
   # LIBS="$LIBS ${netcdf_libs}"

   dnl check if HDF5 version 1.10.4 or later is used to build NetCDF
   if test "x$nc_config" != x ; then
      # cflags should contain HDF5 include path
      netcdf_cflags=`$nc_config --cflags`
   fi

   CPPFLAGS="$CPPFLAGS $netcdf_cflags"
   dnl AC_CHECK_HEADERS uses CPPFLAGS, not CFLAGS
   AC_CHECK_HEADER([hdf5.h], [have_hdf5=yes], [have_hdf5=no])
   if test "x$have_hdf5" = xyes ; then
      AC_MSG_CHECKING([whether HDF5 version is 1.10.4 or later])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <hdf5.h>
#if (H5_VERS_MAJOR*1000000 + H5_VERS_MINOR*1000 + H5_VERS_RELEASE < 1010004)
#error HDF5 version is older than 1.10.4
#endif
      ]])], [hdf5_ge_1_10_4=yes], [hdf5_ge_1_10_4=no])
      AC_MSG_RESULT([$hdf5_ge_1_10_4])
      if test x$hdf5_ge_1_10_4 = xyes; then
         AC_DEFINE([HDF5_VER_GE_1_10_4])
      fi
   fi
   # restore CPPFLAGS
   CPPFLAGS=$saved_CPPFLAGS

   # get NETCDF4 include dir from nc-config
   # AC_MSG_CHECKING([NetCDF-4 include directories])
   # netcdf_incdir=
   # if test "x$nc_config" != x ; then
   #    netcdf_incdir=`$nc_config --includedir`
   # else
   #    if test "x$NETCDF4_INSTALL" != x ; then
   #       netcdf_incdir="$NETCDF4_INSTALL/include"
   #    fi
   # fi
   # NETCDF4_INC=
   # if test "x$netcdf_incdir" = x ; then
   #    AC_MSG_RESULT([none needed])
   # else
   #    AC_MSG_RESULT([$netcdf_incdir])
   #    NETCDF4_INC="-I$netcdf_incdir"
   # fi
   # AC_SUBST([NETCDF4_INC])
   NETCDF4_INC=
   if test "x$netcdf_includedir" != x ; then
      NETCDF4_INC="-I$netcdf_includedir"
   fi
   AC_SUBST([NETCDF4_INC])

   AC_MSG_CHECKING([NetCDF-4 LDFLAGS])
   NETCDF4_LDFLAGS=
   if test "x$netcdf4_libdir" != x ; then
      NETCDF4_LDFLAGS="-L$netcdf4_libdir"
   fi
   UD_MSG_DEBUG(NETCDF4_LDFLAGS=$NETCDF4_LDFLAGS)
   # remove all -l options from netcdf_libs that should leave -L only
   if test "x$netcdf_libs" = x ; then
      netcdf4_L=`echo $netcdf_libs | ${SED} 's/-l[[^ ]]*//g'`
   fi
   if test "x$netcdf4_L" = x ; then
      netcdf4_L=$NETCDF4_LDFLAGS
   elif test "x$NETCDF4_LDFLAGS" != x ; then
      netcdf4_L="$NETCDF4_LDFLAGS $netcdf4_L"
   fi
   # remove space between -L and path
   if test "x$netcdf4_L" != x ; then
      netcdf4_L=`echo $netcdf4_L | ${SED} 's/-L[[ ]]*/-L/g'`
   fi
   UD_MSG_DEBUG(netcdf4_L=$netcdf4_L)

   if test "x$netcdf4_L" != x ; then
      # remove any duplicate in netcdf4_L
      NETCDF4_LDFLAGS=`echo $netcdf4_L | ${AWK} -v RS=' ' '{if (!seen[[$1]]++){printf "%s%s",sep,$1; sep=" "}}'`
   fi
   if test "x$NETCDF4_LDFLAGS" != x ; then
      AC_MSG_RESULT([$NETCDF4_LDFLAGS])
   else
      AC_MSG_RESULT([$LDFLAGS])
   fi
   AC_SUBST([NETCDF4_LDFLAGS])

   # remove all -L options from netcdf_libs that should leave -l only
   netcdf4_l=`echo $netcdf_libs | ${SED} 's/-L[[^ ]]*//g'`
   UD_MSG_DEBUG(netcdf4_l=$netcdf4_l)

   AC_MSG_CHECKING([NetCDF-4 LIBS])
   # remove any duplicate in netcdf4_l
   NETCDF4_LIBS=`echo $netcdf4_l | ${AWK} -v RS=' ' '{if (!seen[[$1]]++){printf "%s%s",sep,$1; sep=" "}}'`
   AC_SUBST([NETCDF4_LIBS])
   AC_MSG_RESULT([$NETCDF4_LIBS])

   AC_SUBST([NETCDF4_PACKAGE], [netcdf])
fi

ENABLE_NETCDF4=0
if test "x$enable_netcdf4" = xyes ; then
   AC_DEFINE(ENABLE_NETCDF4)
   ENABLE_NETCDF4=1
fi
AC_SUBST(ENABLE_NETCDF4)
AM_CONDITIONAL(ENABLE_NETCDF4, [test x$enable_netcdf4 = xyes])
AM_CONDITIONAL(NETCDF_480_481, [test x$netcdf_480_481 = xyes])

AC_ARG_ENABLE([erange-fill],
   [AS_HELP_STRING([--disable-erange-fill],
                   [Disable use of fill value when out-of-range type
                    conversion causes NC_ERANGE error. @<:@default: enabled@:>@])],
   [enable_erange_fill=${enableval}], [enable_erange_fill=yes]
)
ENABLE_ERANGE_FILL=0
if test "x$enable_erange_fill" = xyes ; then
   ENABLE_ERANGE_FILL=1
fi
UD_MSG_DEBUG(enable_erange_fill=$enable_erange_fill)
AC_SUBST(ENABLE_ERANGE_FILL)
AM_CONDITIONAL(ENABLE_ERANGE_FILL, [test x$enable_erange_fill = xyes])

AC_ARG_ENABLE([relax-coord-bound],
   [AS_HELP_STRING([--disable-relax-coord-bound],
                   [Use stricter rule for error NC_EINVALCOORDS to disallow coordinate
                    start argument equal to dimension size when argument
                    count is zero. @<:@default: enabled@:>@])],
   [relax_coord_bound=${enableval}], [relax_coord_bound=auto]
)
UD_MSG_DEBUG(relax_coord_bound=$relax_coord_bound)

if test "x$enable_netcdf4" = xyes ; then
   # check nc_relax_coord_bound set in NetCDF-4
   if test "x$relax_coord_bound" = xauto ; then
      relax_coord_bound=$nc_relax_coord_bound
   else
      if test "$relax_coord_bound" != "$nc_relax_coord_bound" ; then
         # error out when the settings conflict
         if test "x$relax_coord_bound" = xyes ; then
             AC_MSG_ERROR([Enabling relax-coord-bound conflicts with NetCDF-4 setting])
         else
             AC_MSG_ERROR([Disabling relax-coord-bound conflicts with NetCDF-4 setting])
         fi
      fi
      relax_coord_bound=$nc_relax_coord_bound
   fi
else
   if test "x$relax_coord_bound" = xauto ; then
      # default is yes
      relax_coord_bound=yes
   fi
fi
UD_MSG_DEBUG(relax_coord_bound=$relax_coord_bound)

RELAX_COORD_BOUND=0
if test "x$relax_coord_bound" = xyes ; then
   RELAX_COORD_BOUND=1
   AC_DEFINE(RELAX_COORD_BOUND)
fi
AC_SUBST(RELAX_COORD_BOUND)
AM_CONDITIONAL([RELAX_COORD_BOUND], [test "x$relax_coord_bound" = xyes])

AC_ARG_ENABLE([doxygen],
   [AS_HELP_STRING([--enable-doxygen], [Enable generation of documentation. @<:@default: disabled@:>@])],
   [enable_doxygen=${enableval}], [enable_doxygen=no]
)
AM_CONDITIONAL([BUILD_DOCS], [test "x$enable_doxygen" = xyes])

AC_PATH_PROG([LATEX],  [latex])
AC_PATH_PROG([DVIPDF], [dvipdf])
has_latex=no
if test "x${LATEX}" != x ; then
    has_latex=yes
fi
AC_SUBST(LATEX)
AC_SUBST(DVIPDF)
AM_CONDITIONAL([HAS_LATEX], [test "x$has_latex" = xyes])

AC_ARG_ENABLE([file-sync],
    [AS_HELP_STRING([--disable-file-sync],
                    [Disable MPI file sync if you know your file system can
                     provide data consistency. @<:@default: enabled@:>@])],
    [file_sync=${enableval}], [file_sync=yes]
)
if test "x${file_sync}" = xno ; then
    AC_DEFINE(DISABLE_FILE_SYNC)
fi

AC_ARG_ENABLE([large-single-req],
    [AS_HELP_STRING([--enable-large-single-req],
                    [Enable large (> 2 GiB) single request in individual MPI-IO
                     calls. Note some MPI-IO libraries may not support this.
                     @<:@default: disabled@:>@])],
    [large_single_req=${enableval}], [large_single_req=no]
)
ENABLE_LARGE_SINGLE_REQ=0
if test "x${large_single_req}" = xyes ; then
    ENABLE_LARGE_SINGLE_REQ=1
    AC_DEFINE(ENABLE_LARGE_SINGLE_REQ)
fi
AM_CONDITIONAL(ENABLE_LARGE_SINGLE_REQ, [test x$large_single_req = xyes])
AC_SUBST(ENABLE_LARGE_SINGLE_REQ)

AC_ARG_ENABLE([large-file-test],
    [AS_HELP_STRING([--enable-large-file-test],
                    [Enable testing for large (>4GB) file/variable I/O. Note
                     "make check" can run very slow. @<:@default: disabled@:>@])],
    [large_file_test=${enableval}], [large_file_test=no]
)
AM_CONDITIONAL(RUN_LARGE_FILE_TEST, [test x$large_file_test = xyes])

# Check whether to enable strict null byte header padding.
# See https://github.com/Unidata/netcdf-c/issues/657 for more information.
AC_ARG_ENABLE([null-byte-header-padding],
    [AS_HELP_STRING([--enable-null-byte-header-padding],
                    [Enable check for null-byte header padding when reading
                     files in classic formats. @<:@default: disabled@:>@])],
    [null_byte_header_padding=${enableval}], [null_byte_header_padding=no]
)
ENABLE_NULL_BYTE_HEADER_PADDING=0
if test "x${null_byte_header_padding}" = xyes ; then
    ENABLE_NULL_BYTE_HEADER_PADDING=1
    AC_DEFINE(ENABLE_NULL_BYTE_HEADER_PADDING)
fi
AC_SUBST(ENABLE_NULL_BYTE_HEADER_PADDING)

AC_ARG_ENABLE([burst-buffering],
   [AS_HELP_STRING([--enable-burst-buffering],
                   [Enable burst buffer driver support. @<:@default: disabled@:>@])],
   [enable_bbdriver=${enableval}], [enable_bbdriver=no]
)

ENABLE_BURST_BUFFER=0
if test "x$enable_bbdriver" = "xyes" ; then
   AC_DEFINE(ENABLE_BURST_BUFFER)
   ENABLE_BURST_BUFFER=1
fi
AC_SUBST(ENABLE_BURST_BUFFER)
AM_CONDITIONAL(ENABLE_BURST_BUFFER, [test x$enable_bbdriver = xyes])

ADIOS_INSTALL=""
AC_ARG_WITH(adios,
   [AS_HELP_STRING([--with-adios@<:@=DIR@:>@],
                   [Enable ADIOS feature and provide the ADIOS installation path in DIR.
                    @<:@default: disabled@:>@])],
   [ dnl this clause is run when --with-adios or --without-adios is used
   if test "x${withval}" = xyes ; then
      enable_adios=yes
   elif test "x${withval}" != xno ; then
      ADIOS_INSTALL=${withval}
      enable_adios=yes
   fi
   ],[enable_adios=no]
)

ENABLE_ADIOS=0
if test "x$enable_adios" = "xyes" ; then
   AC_DEFINE(ENABLE_ADIOS)
   ENABLE_ADIOS=1
fi
AC_SUBST(ENABLE_ADIOS)
AM_CONDITIONAL(ENABLE_ADIOS, [test x$enable_adios = xyes])

ADIOS_VER_GE_1132=0
if test "x$enable_adios" = "xyes" ; then
   if test "x$ADIOS_INSTALL" = x ; then
      AC_PATH_PROG([adios_config],[adios_config],,[$PATH])
   else
      AC_PATH_PROG([adios_config],[adios_config],,[$ADIOS_INSTALL/bin])
   fi

   if test "x$adios_config" = x ; then
      AC_MSG_ERROR([
      ------------------------------------------------------------
      Cannot find adios_config, a utility command to show the
      configuration of the ADIOS library. It is required to
      build PnetCDF with BP read support. Use option
          --with-adios=/path/to/implementation
      to specify the location of ADIOS build. In addition,
      please make sure the MPI C compiler is compatible with the
      one used to build ADIOS. The ADIOS library used must
      be built with parallel I/O support. Check
      'config.log' for more information. Stopping ...
      ------------------------------------------------------------])
   fi
   UD_MSG_DEBUG([adios_config=$adios_config])

   # Get adios compile flags
   adios_cflag=`$adios_config -c`
   UD_MSG_DEBUG([adios_compiler_flag=$adios_cflag])

   # Get adios linker flags
   adios_ldflag=`$adios_config -l`
   UD_MSG_DEBUG([adios_linker_flag=$adios_ldflag])

   # remove all -l options from adios_ldflag that should leave -L only
   adios_LDFLAGS=`echo $adios_ldflag | ${SED} 's/-l[[^ ]]*//g'`
   UD_MSG_DEBUG(after removing -l adios_LDFLAGS=$adios_LDFLAGS)

   # remove space between -L and path
   adios_LDFLAGS=`echo $adios_LDFLAGS | ${SED} 's/-L[[ ]]*/-L/g'`
   UD_MSG_DEBUG(after removing space between -L and path adios_LDFLAGS=$adios_LDFLAGS)

   # remove -L from adios_LDFLAGS that should leave libdir only
   adios_libdir=`echo ${adios_LDFLAGS} | ${SED} 's/-L//g'`
   UD_MSG_DEBUG(adios_libdir=$adios_libdir)

   # remove all -L options from adios_ldflag that should leave -l only
   adios_LIBS=`echo $adios_ldflag | ${SED} 's/-L[[^ ]]*//g'`
   UD_MSG_DEBUG(adios_LIBS=$adios_LIBS)

   saved_CPPFLAGS=$CPPFLAGS
   saved_LDFLAGS=$LDFLAGS
   saved_LIBS=$LIBS

   CPPFLAGS="$adios_cflag $CPPFLAGS"
   LDFLAGS="$adios_LDFLAGS $LDFLAGS"
   LIBS="$adios_LIBS $LIBS"

   have_adios=no
   AC_SEARCH_LIBS([adios_read_open_file], [adiosread], [have_adios=yes], [have_adios=no])
   if test "x${have_adios}" = xyes; then
      AC_SEARCH_LIBS([adios_posix_open_read_internal], [adios], [have_adios=yes], [have_adios=no])
      if test "x${have_adios}" = xyes; then
         AC_CHECK_HEADERS([adios_read.h], [], [have_adios=no])
      fi
   fi

   if test "x${have_adios}" = xno; then
      have_adios2=no
      AC_MSG_CHECKING(ADIOS2 library)
      AC_CHECK_HEADERS([adios2_c.h], [have_adios2=yes], [have_adios2=no])
      if test "x${have_adios2}" = xno; then
         AC_MSG_ERROR([
         ------------------------------------------------------------
         The ADIOS library and header file are required to build
         PnetCDF with ADIOS BP read support. Use option
             --with-adios=/path/to/implementation
         to specify the location of ADIOS build. In addition,
         please make sure the MPI C compiler is compatible with the
         one used to build ADIOS.
         Stopping ...
         Check 'config.log' for more information.
         ------------------------------------------------------------])
      else
         AC_MSG_ERROR([
         ------------------------------------------------------------
         PnetCDF currently only supports ADIOS 1.x.
         ADIOS 1.x can be downloaded at:
         https://github.com/ornladios/ADIOS
         We will support ADIOS 2.x in future releases.
         ------------------------------------------------------------])
      fi
   else
      LIBS="$LIBS -ladiosread -ladios"
   fi

   dnl Check if ADIOS version is 1.12.0 or later
   AC_MSG_CHECKING([Check if ADIOS version is 1.12.0 or later])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <adios_version.h>
#if (ADIOS_VERSION_MAJOR*1000000 + ADIOS_VERSION_MINOR*1000 + ADIOS_VERSION_PATCH < 1012000)
      choke me
#endif
   ]])], [adios_ge_1_12_0=yes], [adios_ge_1_12_0=no])
   AC_MSG_RESULT([$adios_ge_1_12_0])

   if test "x${adios_ge_1_12_0}" = xno; then
      AC_MSG_ERROR([
      ------------------------------------------------------------
      ADIOS version 1.12.0 or later is required to build
      PnetCDF with ADIOS BP read support.
      ------------------------------------------------------------])
   fi

   dnl Check if ADIOS version is 1.13.2 or later
   AC_MSG_CHECKING([Check if ADIOS version is 1.13.2 or later])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <adios_version.h>
#if (ADIOS_VERSION_MAJOR*1000000 + ADIOS_VERSION_MINOR*1000 + ADIOS_VERSION_PATCH < 1013002)
      choke me
#endif
   ]])], [adios_ge_1_13_2=yes], [adios_ge_1_13_2=no])
   AC_MSG_RESULT([$adios_ge_1_13_2])

   if test "x${adios_ge_1_13_2}" = xyes; then
      ADIOS_VER_GE_1132=1
   fi

   CPPFLAGS=$saved_CPPFLAGS
   LDFLAGS=$saved_LDFLAGS
   LIBS=$saved_LIBS

   AC_SUBST([ADIOS_INC], ["$adios_cflag"])
   AC_SUBST([ADIOS_LDFLAGS], ["$adios_LDFLAGS"])
   AC_SUBST([ADIOS_LIBS], ["$adios_LIBS"])
   AC_SUBST([ADIOS_PACKAGE], [adios])
fi

AC_SUBST(ADIOS_VER_GE_1132)

AC_ARG_ENABLE([profiling],
   [AS_HELP_STRING([--enable-profiling],
                   [Enable time and memory profiling. @<:@default: disabled@:>@])],
   [enable_profiling=${enableval}], [enable_profiling=no]
)

PNETCDF_PROFILING=0
if test "x$enable_profiling" = "xyes" ; then
   AC_DEFINE(PNETCDF_PROFILING)
   PNETCDF_PROFILING=1

   if test "x${debug}" = xno; then
      dnl check required functions for enabling malloc tracing
      if test "x${ac_cv_func_tsearch}" = xyes && test "x${ac_cv_func_tdelete}" = xyes ; then
         AC_DEFINE(PNC_MALLOC_TRACE)
      fi
   fi
fi
AC_SUBST(PNETCDF_PROFILING)
AM_CONDITIONAL(PNETCDF_PROFILING, [test x$enable_profiling = xyes])

dnl build test programs and benchmark programs
AM_CONDITIONAL(BUILD_TESTSETS, [true])
AM_CONDITIONAL(BUILD_BENCHMARKS_IN_PNETCDF, [true])

AC_ARG_VAR(TESTSEQRUN, [Run command (on one MPI process) for "make check" on cross-compile environment. Example: "aprun -n 1". @<:@default: none@:>@])

AC_ARG_VAR(TESTMPIRUN, [MPI run command for "make ptest", @<:@default: mpiexec -n NP@:>@])
if test "x${TESTMPIRUN}" = x ; then
   dnl if TESTMPIRUN has not been set by users, then
   dnl set default to "mpiexec -n NP"
   UD_MPI_PATH_PROGS([TESTMPIRUN], [mpiexec mpirun srun])
   if test "x${TESTMPIRUN}" != x ; then
      TESTMPIRUN="$TESTMPIRUN -n NP"
   fi
fi

AC_ARG_VAR(TESTOUTDIR, [Output file directory for "make check" and "make ptest", @<:@default: ./@:>@])
if test "x${TESTOUTDIR}" = x ; then
   dnl set default to current directory
   TESTOUTDIR=.
fi

# SEQ_CC is used to compile programs to be run sequentially, such as
# pnetcdf_version, ncoffsets, and ncvalidator
# Note on Theta @ALCF, CPATH is set for Intel Compiler and when gcc is used to
# compile the above utility programs, errors can occur. See doc/README.CRAY for
# more information. Below we use the same C compiler as MPICC.
AC_ARG_VAR(SEQ_CC, [C compiler to build sequential utility programs, @<:@default: gcc@:>@])
if test "x$SEQ_CC" = x ; then
   # if test "x$ac_cv_mpicc_base" = xGCC ; then
   #    AC_PATH_PROG([SEQ_CC], [gcc], [$MPICC])
   # elif test "x$ac_cv_mpicc_base" = xCLANG ; then
   #    AC_PATH_PROG([SEQ_CC], [clang], [$MPICC])
   # elif test "x$ac_cv_mpicc_base" = xICC ; then
   #    AC_PATH_PROG([SEQ_CC], [icc], [$MPICC])
   # elif test "x$ac_cv_mpicc_base" = xPGCC ; then
   #    AC_PATH_PROG([SEQ_CC], [pgcc], [$MPICC])
   # elif test "x$ac_cv_mpicc_base" = xXLC ; then
   #    AC_PATH_PROG([SEQ_CC], [xlc], [$MPICC])
   # elif test "x$ac_cv_mpicc_base" = xFCCPX ; then
   #    AC_PATH_PROG([SEQ_CC], [fccpx], [$MPICC])
   # fi
   if test "x$ac_cv_mpi_compiler_base_MPICC" = x ; then
      # gcc should be available on most of the machines today
      AC_PATH_PROG([SEQ_CC], [gcc], [$MPICC])
   else
      SEQ_CC=$ac_cv_mpi_compiler_base_MPICC
   fi
fi
AC_MSG_CHECKING([C compiler for serial utility programs])
AC_MSG_RESULT([$SEQ_CC])
AC_SUBST(SEQ_CC)

dnl Configuration Date
dnl Note that command 'date' is not portable across Unix platforms
if test "x$SOURCE_DATE_EPOCH" != x ; then
    date_str=`date -d @${SOURCE_DATE_EPOCH}`
    if test "x$date_str" != x ; then
       AC_SUBST([CONFIG_DATE], ["$date_str"])
    else
       date_str=`date -r $SOURCE_DATE_EPOCH`
       if test "x$date_str" != x ; then
          AC_SUBST([CONFIG_DATE], ["$date_str"])
       else dnl ignore SOURCE_DATE_EPOCH
          AC_SUBST([CONFIG_DATE], ["`date`"])
       fi
    fi
else
    AC_SUBST([CONFIG_DATE], ["`date`"])
fi
UD_MSG_DEBUG([CONFIG_DATE=$CONFIG_DATE])

dnl Enable creation of libtool-style versioning or no versioning
AC_ARG_ENABLE(versioning,
        [AC_HELP_STRING([--disable-versioning],[Disable library versioning. @<:@default: enabled@:>@])],
        [enable_versioning=${enableval}], [enable_versioning=yes]
)
dnl For libtool ABI versioning rules see:
dnl http://www.gnu.org/software/libtool/manual/libtool.html#Updating-version-info
dnl Update the version information only immediately before a public release.
dnl PnetCDF starts with 1:0:0 (shared library is first supported in 1.9.0)
dnl because some package distributors, such as Debian, may have already built
dnl PnetCDF with shared libraries.
ABIVERSION="4:3:0"
AC_SUBST(ABIVERSION)
if test "$enable_versioning" = "yes" ; then
   ABIVERSIONFLAGS="-version-info \$(ABIVERSION)"
else
   ABIVERSIONFLAGS="-avoid-version"
fi
AC_SUBST(ABIVERSIONFLAGS)

dnl AC_CONFIG_HEADERS([src/binding/f77/nfconfig_inc])
AC_CONFIG_FILES(Makefile \
                src/Makefile \
                src/dispatchers/Makefile \
                src/libs/Makefile \
                src/include/Makefile \
                src/include/pnetcdf.h \
                src/drivers/Makefile \
                src/drivers/common/Makefile \
                src/drivers/include/Makefile \
                src/drivers/ncmpio/Makefile \
                src/drivers/nc4io/Makefile \
                src/drivers/ncadios/Makefile \
                src/drivers/ncbbio/Makefile \
                src/drivers/ncfoo/Makefile \
                src/binding/Makefile \
                src/binding/cxx/Makefile \
                src/binding/f77/Makefile \
                src/binding/f77/pnetcdf.inc \
                src/binding/f90/Makefile \
                src/binding/f90/pnetcdf.f90 \
                src/binding/f90/api.fh \
                src/binding/f90/nfmpi_constants.fh \
                src/utils/Makefile \
                src/utils/ncoffsets/Makefile \
                src/utils/pnetcdf_version/Makefile \
                src/utils/ncvalidator/Makefile \
                src/utils/ncmpidiff/Makefile \
                src/utils/ncmpidump/Makefile \
                src/utils/ncmpigen/Makefile \
                src/utils/ncmpilogdump/Makefile \
                src/utils/pnetcdf-config \
                src/packaging/Makefile \
                src/packaging/pnetcdf.pc \
                examples/Makefile \
                examples/adios/Makefile \
                examples/C/Makefile \
                examples/tutorial/Makefile \
                examples/CXX/Makefile \
                examples/F77/Makefile \
                examples/F90/Makefile \
                examples/burst_buffer/Makefile \
                benchmarks/Makefile \
                benchmarks/C/Makefile \
                benchmarks/FLASH-IO/Makefile \
                doc/Makefile \
                doc/pnetcdf-api/Makefile \
                man/Makefile \
                test/Makefile \
                test/common/Makefile \
                test/C/Makefile \
                test/fandc/Makefile \
                test/nc_test/Makefile \
                test/cdf_format/Makefile \
                test/nc4/Makefile \
                test/adios/Makefile \
                test/header/Makefile \
                test/testcases/Makefile \
                test/nonblocking/Makefile \
                test/largefile/Makefile \
                test/CXX/Makefile \
                test/subfile/Makefile \
                test/nf_test/Makefile \
                test/nf_test/tests.inc \
                test/nf90_test/Makefile \
                test/nf90_test/tests.inc \
                test/F90/Makefile \
                test/burst_buffer/Makefile \
)

dnl The following dependency is for configure.in and configure
dnl See autoconf manual 2.69, Section 4.8.5 Automatic Remaking
dnl AC_CONFIG_FILES([stamp-h], [echo timestamp > stamp-h])

dnl build ncfoo driver only in debug mode
AM_CONDITIONAL(BUILD_DRIVER_FOO, [test x${debug} = xyes])

AC_SUBST(enable_static)
AC_SUBST(enable_shared)
AC_SUBST(enable_netcdf4)
AC_SUBST(netcdf4_libdir)
AC_SUBST(enable_adios)
AC_SUBST(adios_libdir)
AM_CONDITIONAL(BUILD_SHARED_LIB, [test x${enable_shared} = xyes])

AC_OUTPUT
echo "------------------------------------------------------------------------------"
echo \
"
   ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

   Features:  Build static libraries                      - ${enable_static}
              Build shared libraries                      - ${enable_shared}
              Build Fortran APIs                          - ${has_fortran}
              Build C++ APIs                              - ${has_mpicxx}"
if test "x${enable_netcdf4}" = xyes; then
   echo "\
              NetCDF-4 support                            - enabled"
fi
if test "x${enable_adios}" = xyes; then
   echo "\
              ADIOS BP read support                       - enabled"
fi
if test "x${enable_bbdriver}" = xyes ; then
   echo "\
              Burst buffering                             - enabled"
fi
if test "x${enable_subfiling}" = xyes; then
   echo "\
              Subfiling support                           - enabled"
fi
if test "x${large_single_req}" = xyes ; then
   echo "\
              Allow large (> 2GiB) single MPI-IO requests - enabled"
fi
if test "x${thread_safe}" = xyes ; then
   echo "\
              Thread-safe capability                      - enabled"
fi
if test "x${enable_aggreg}" = xno; then
   echo "\
              Request aggregation in nonblocking APIs     - disabled"
fi
if test "x${enable_erange_fill}" = xno; then
   echo "\
              Fill variables when NC_ERANGE occurs        - disabled"
fi
if test "x${relax_coord_bound}" = xno; then
   echo "\
              Relax start coordinate bound check          - disabled"
fi
if test "x${ac_cv_c_bigendian}" = xno ; then
   if test "x${in_place_swap}" = xno ; then
      echo "\
              Memory in-place byte swap                   - disabled"
   elif test "x${in_place_swap}" = xyes ; then
      echo "\
              Memory in-place byte swap                   - enabled"
   fi
fi
if test "x${large_file_test}" = xyes; then
   echo "\
              Testing large file/variable I/O             - enabled"
fi
if test "x${has_fortran}" = xyes && test "x${ac_f77_support_freeform}" = xno ; then
   echo "\
              Support free form in Fortran 77             - disabled"
fi
if test "x${null_byte_header_padding}" = xyes ; then
   echo "\
              Check null-byte header padding for reads    - enabled"
fi
if test "x${debug}" = xyes; then
   echo "\
              Internal debug mode                         - enabled"
fi
if test "x${enable_profiling}" = xyes ; then
   echo "\
              Internal profiling mode                     - enabled"
fi

echo "\

   Compilers: MPICC    = ${MPICC}"
if test "${has_mpicxx}" = yes ; then
   echo "\
              MPICXX   = ${MPICXX}"
fi
if test "${has_fortran}" = yes ; then
   echo "\
              MPIF77   = ${MPIF77}
              MPIF90   = ${MPIF90}"
fi
if test "x${CPPFLAGS}" != x ; then
   echo "\
              CPPFLAGS = ${CPPFLAGS}"
fi
echo "\
              CFLAGS   = ${CFLAGS}"
if test "${has_mpicxx}" = yes ; then
   echo "\
              CXXFLAGS = ${CXXFLAGS}"
fi
if test "${has_fortran}" = yes ; then
   echo "\
              FFLAGS   = ${FFLAGS}
              FCFLAGS  = ${FCFLAGS}"
fi
if test "x${LDFLAGS}" != x ; then
   echo "\
              LDFLAGS  = ${LDFLAGS}"
fi
if test "x${LIBS}" != x ; then
   echo "\
              LIBS     = ${LIBS}"
fi
if test "x${enable_netcdf4}" = xyes; then
   echo "\

   NetCDF-4 library:
              NETCDF4_INC     = ${NETCDF4_INC}
              NETCDF4_LDFLAGS = ${NETCDF4_LDFLAGS}
              NETCDF4_LIBS    = ${NETCDF4_LIBS}"
fi
if test "x${enable_adios}" = xyes; then
   echo "\

   ADIOS library:
              ADIOS_INC     = ${ADIOS_INC}
              ADIOS_LDFLAGS = ${ADIOS_LDFLAGS}
              ADIOS_LIBS    = ${ADIOS_LIBS}"
fi
echo "\

   Now run 'make' to build the library and utility tools.
   Then run 'make @<:@<target>@:>@' for testing and installation, where the
   optional <target> can be:
              tests    - build all test programs (build only, no run)
              check    - run sequential test programs
              ptest    - run parallel test programs on 4 MPI processes
              ptests   - run parallel test programs on 3,4,6,8 MPI processes
              install  - install PnetCDF library in ${prefix}

------------------------------------------------------------------------------"

