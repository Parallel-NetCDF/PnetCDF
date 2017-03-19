AC_REVISION($Id$)dnl
dnl -*- Mode: shell-script-mode; -*-
dnl Process this file with GNU autoconf(1) to produce a configure script.
dnl

dnl autoconf v2.69 was released in 2012-04-24
AC_PREREQ([2.69])
AC_INIT([parallel-netcdf],[1.8.0],[parallel-netcdf@mcs.anl.gov])

dnl ncconfig.h.in will be created by autoreconf (autoheader)
dnl call it right after AC_INIT, as suggested by autoconf
AC_CONFIG_HEADERS([src/include/ncconfig.h])

AC_CONFIG_SRCDIR([src/include/pnetcdf.h.in])

AC_CONFIG_AUX_DIR([./scripts])

dnl Note getting command line should be done before calling AM_INIT_AUTOMAKE
dnl as AM_INIT_AUTOMAKE modfies command line $*
CONFIGURE_ARGS_CLEAN=`echo $* | tr '"' ' '`

dnl AM_INIT_AUTOMAKE([subdir-objects])
AM_INIT_AUTOMAKE([1.13])
dnl enable silent rules by default 
AM_SILENT_RULES([yes])

dnl if maintainer mode is disabled, make will *never* attempt to rebuild configure,
dnl Makefile.ins, etc.
AM_MAINTAINER_MODE([enable])

m4_ifdef([AC_CONFIG_MACRO_DIRS], [AC_CONFIG_MACRO_DIRS([m4])], [AC_CONFIG_MACRO_DIR([m4])])

dnl AM_EXTRA_RECURSIVE_TARGETS macro was introduced into automake 1.13
m4_ifdef([AM_EXTRA_RECURSIVE_TARGETS], [AM_EXTRA_RECURSIVE_TARGETS([tests])])

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
SVN_DATE="$LastChangedDate$"
PNETCDF_RELEASE_DATE2=`echo $SVN_DATE | cut -d' ' -f2`
PNETCDF_RELEASE_DATE=`echo $SVN_DATE | cut -d' ' -f6,7,8 | cut -d')' -f1`

dnl user defined macro for printing messages for debugging
_DEBUG=no
AC_DEFUN([UD_MSG_DEBUG],
         [if test "x${_DEBUG}" = xyes ; then
             AC_MSG_NOTICE([DEBUG: $1])
          fi
         ]
)

UD_MSG_DEBUG([PNETCDF_VERSION_MAJOR=$PNETCDF_VERSION_MAJOR])
UD_MSG_DEBUG([PNETCDF_VERSION_MINOR=$PNETCDF_VERSION_MINOR])
UD_MSG_DEBUG([PNETCDF_VERSION_SUB=$PNETCDF_VERSION_SUB])
UD_MSG_DEBUG([PNETCDF_VERSION_PRE=$PNETCDF_VERSION_PRE])
UD_MSG_DEBUG([PNETCDF_VERSION=$PNETCDF_VERSION])
UD_MSG_DEBUG([PNETCDF_RELEASE_DATE=$PNETCDF_RELEASE_DATE])

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
AC_SUBST(PNETCDF_RELEASE_DATE)
AC_SUBST(PNETCDF_RELEASE_DATE2)
AC_SUBST(PACKAGE_VERSION)
AC_SUBST(CONFIGURE_ARGS_CLEAN)

dnl these used to live in acconfig.h
dnl autoheader only adds these templates to the first invocation of AC_CONFIG_HEADERS
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
AH_TEMPLATE([DISABLE_IN_PLACE_SWAP],    [Define if to disable in-place byte swap])
AH_TEMPLATE([ENABLE_SUBFILING],         [Define if to enable subfiling feature])
AH_TEMPLATE([PNC_MALLOC_TRACE],         [Define if to enable malloc tracing])
AH_TEMPLATE([HAVE_FUNC_MACRO],          [Define if C++ macro __func__ is defined])
AH_TEMPLATE([HAVE_FUNCTION_MACRO],      [Define if C++ macro __FUNCTION__ is defined])
AH_TEMPLATE([RELAX_COORD_BOUND],        [Define if relaxed coord check is enabled])

dnl an option to use a customized rm command
AC_ARG_VAR(RM, Command for deleting files or directories. @<:@default: rm@:>@)
if test "x${RM}" != x ; then
   AC_MSG_CHECKING(rm )
   if ! test -f ${RM} ; then
      AC_CHECK_PROG([rm_cmd], [${RM}], [yes], [no])
      if test "x${rm_cmd}" = xyes ; then
         RM=${RM}
      fi
   else
      RM=${RM}
   fi
   AC_MSG_RESULT(using $RM)
else
   RM="rm"
fi
AC_SUBST(RM)

AC_ARG_ENABLE(echo,
   [AS_HELP_STRING([--enable-echo],
                   [Turn on strong echoing. @<:@default: disabled@:>@])],
   [set -x]
)

MPI_INSTALL=
AC_ARG_WITH(mpi,
   [AS_HELP_STRING([--with-mpi=/path/to/implementation],
                   [installation prefix for MPI implementation])],
   if test x"$withval" = xyes; then
      AC_MSG_ERROR(--with-mpi must be given a pathname)
   else
      MPI_INSTALL=${withval}
   fi
)
if test "x${MPI_INSTALL}" != x && (! test -d "${MPI_INSTALL}") ; then
   AC_MSG_ERROR(Directory '${MPI_INSTALL}' specified in --with-mpi does not exist)
fi

AC_ARG_VAR(MPICC,  MPI C compiler)
AC_ARG_VAR(MPICXX, MPI C++ compiler)
AC_ARG_VAR(MPIF77, MPI Fortran 77 compiler)
AC_ARG_VAR(MPIF90, MPI Fortran 90 compiler)
AC_ARG_VAR(CC,     Overwritten by MPICC  if MPICC is set)
AC_ARG_VAR(CXX,    Overwritten by MPICXX if MPICXX is set)
AC_ARG_VAR(F77,    Overwritten by MPIF77 if MPIF77 is set)
AC_ARG_VAR(FC,     Overwritten by MPIF90 if MPIF90 is set)
dnl AC_ARG_VAR(CPPFLAGS, [Preprocessor options for C and C++ compilers, e.g. -I<include_dir> if you have headers in a nonstandard directory <include_dir>])
dnl AC_ARG_VAR(CFLAGS, Debugging and optimization options for the C compiler)
dnl AC_ARG_VAR(CXXFLAGS, Debugging and optimization options for the C++ compiler)
dnl AC_ARG_VAR(FFLAGS, Debugging and optimization options for the Fortran 77 compiler)
dnl AC_ARG_VAR(FCFLAGS, Debugging and optimization options for the Fortran 90 compiler)

dnl Check if MPICC, MPICXX, MPIF77, MPIF90 are set by the user.
dnl If not, set MPICC equal to CC. Similarly, for MPICXX, MPIF77, and MPIF90
if test "x$MPICC"  = x && test "x$CC"  != x ; then MPICC=$CC   ; fi
if test "x$MPICXX" = x && test "x$CXX" != x ; then MPICXX=$CXX ; fi
if test "x$MPIF77" = x && test "x$F77" != x ; then MPIF77=$F77 ; fi
if test "x$MPIF77" = x && test "x$FC"  != x ; then MPIF77=$FC  ; fi
if test "x$MPIF90" = x && test "x$F90" != x ; then MPIF90=$F90 ; fi
if test "x$MPIF90" = x && test "x$FC"  != x ; then MPIF90=$FC  ; fi

CANDIDATE_MPICC="${MPICC} mpicc"
CANDIDATE_MPICXX="${MPICXX} mpicxx mpic++ mpiCC"
CANDIDATE_MPIF77="${MPIF77} mpif77"
CANDIDATE_MPIF90="${MPIF90} mpif90"

dnl add IBM MPI compilers
CANDIDATE_MPICC="${CANDIDATE_MPICC} mpcc_r mpcc mpixlc_r mpixlc"
CANDIDATE_MPICXX="${CANDIDATE_MPICXX} mpCC_r mpCC mpixlcxx_r mpixlcxx mpixlC_r mpixlC"
CANDIDATE_MPIF77="${CANDIDATE_MPIF77} mpixlf77_r mpixlf77"
CANDIDATE_MPIF90="${CANDIDATE_MPIF90} mpixlf90_r mpixlf90"

dnl add IBM BGL MPI compilers
CANDIDATE_MPICC="${CANDIDATE_MPICC} blrts_xlc mpxlc_r mpxlc"
CANDIDATE_MPICXX="${CANDIDATE_MPICXX} blrts_xlC mpxlC_r mpxlC"
CANDIDATE_MPIF77="${CANDIDATE_MPIF77} blrts_xlf mpxlf_r mpxlf"
CANDIDATE_MPIF90="${CANDIDATE_MPIF90} blrts_xlf90 mpxlf90_r mpxlf90 mpxlf95_r mpxlf95"

dnl add Fujitsu MPI compilers
CANDIDATE_MPICC="${CANDIDATE_MPICC} mpifccpx"
CANDIDATE_MPICXX="${CANDIDATE_MPICXX} mpiFCCpx"
CANDIDATE_MPIF77="${CANDIDATE_MPIF77} mpifrtpx"
CANDIDATE_MPIF90="${CANDIDATE_MPIF90} mpifrtpx"

dnl add Cray MPI compiler wrappers
CANDIDATE_MPICC="${CANDIDATE_MPICC} cc"
CANDIDATE_MPICXX="${CANDIDATE_MPICXX} CC"
CANDIDATE_MPIF77="${CANDIDATE_MPIF77} ftn"
CANDIDATE_MPIF90="${CANDIDATE_MPIF90} ftn"

dnl add Intel MPI compiler wrappers
CANDIDATE_MPICC="${CANDIDATE_MPICC} mpiicc icc"
CANDIDATE_MPICXX="${CANDIDATE_MPICXX} mpiicpc icpc"
CANDIDATE_MPIF77="${CANDIDATE_MPIF77} mpiifort mpiifc ifort"
CANDIDATE_MPIF90="${CANDIDATE_MPIF90} mpiifort mpiifc ifort"

dnl if MPICC or CC has not been set by users, then search from CANDIDATE_MPICC
if test "x${MPICC}" = x ; then
   UD_MPI_PATH_PROGS([MPICC], [$CANDIDATE_MPICC])
fi
if test "x${MPICC}" = x ; then
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
    No MPI C compiler can be found. Parallel netCDF requires an MPI C
    compiler. Please specify the location of one either with the MPICC
    environment variable or the --with-mpi configure flag
   -----------------------------------------------------------------------])
fi

CC=${MPICC}
AC_SUBST(MPICC)

dnl Now MPICC is configured, do some basic compiler tests
AC_PROG_CC([${MPICC}])
AM_PROG_CC_C_O

dnl this call needs at least autoconf version 2.60
dnl AC_USE_SYSTEM_EXTENSIONS

dnl enable large file support
AC_SYS_LARGEFILE

dnl UD_PROG_CC_MAKEDEPEND

dnl AC_HEADER_STDC
dnl AC_CHECK_HEADERS([malloc.h])

# check if MPICC works for basic MPI call: MPI_Comm_rank()
AC_CHECK_FUNC([MPI_Comm_rank], [],
   dnl maybe -lmpi is needed at link stage
   [AC_SEARCH_LIBS([MPI_Comm_rank], [mpi mpich], [],
                 [AC_MSG_ERROR([
    ------------------------------------------------------------
     Invalid MPI compiler specified or detected: "${MPICC}"
     A working MPI compiler is required. Please specify the location
     of one either with the MPICC environment variable or the
     --with-mpi configure flag
    ------------------------------------------------------------])
])])


AC_ARG_ENABLE([mpi-io-test],
   [AS_HELP_STRING([--disable-mpi-io-test],
                   [Disable check for MPI-IO support in MPI implementation,
                    if you know your MPI implementation has MPI-IO support
                    but the configure test fails to find it.
                    @<:@default: enabled@:>@])],
   [enable_mpi_io_test=${enableval}], [enable_mpi_io_test=yes]
)

if test "x${enable_mpi_io_test}" = xyes ; then
   AC_CHECK_FUNC([MPI_File_open], [],
      dnl maybe -lmpio is needed at link stage
      [AC_SEARCH_LIBS([MPI_File_open], [mpio], [],
                    [AC_MSG_ERROR([
            ------------------------------------------------------------
             MPI implementation does not support MPI-IO
             PnetCDF requires MPI-IO support to work properly.
            ------------------------------------------------------------])]
   )])
fi

dnl check MPI C++ compiler
AC_ARG_ENABLE(cxx,
    [AS_HELP_STRING([--disable-cxx],
                    [Turn off support for the C++ interface,
                     if you only need the C interface. @<:@default: enabled@:>@])],
    [enable_cxx=${enableval}], [enable_cxx=auto]
)
UD_MSG_DEBUG(enable_cxx=$enable_cxx)

if test "x${enable_cxx}" != xno ; then
   dnl if MPICXX or CXX has not been set by users, then search from
   dnl CANDIDATE_MPICXX
   if test "x${MPICXX}" = x ; then
      UD_MPI_PATH_PROGS([MPICXX], [$CANDIDATE_MPICXX])
   fi

   if test "x${MPICXX}" = x ; then
      AC_MSG_ERROR([
         ------------------------------------------------------------
           configure cannot find a C++ compiler. Please specify the
           locations of the compiler either with the MPICXX
           environment variable or the --with-mpi configure flag.
         ------------------------------------------------------------])
      has_mpicxx=no
   else
      CXX=${MPICXX}
      AC_SUBST(MPICXX)
      has_mpicxx=yes

      dnl UD_PROG_CXX
      AC_PROG_CXX(${MPICXX})
      dnl autoconf 2.5.9 has not yet implemented AC_PROG_CXX_C_O
      dnl AC_PROG_CXX_C_O

      dnl test if MPICXX can compile an MPI-IO program
      if test "x${enable_mpi_io_test}" = xyes ; then
          AC_LANG_PUSH(C++)
          AC_CHECK_FUNC([MPI_File_close], [],
             dnl maybe -lmpi++ is needed at link stage
             [AC_SEARCH_LIBS([MPI_File_close], [mpi++ mpichcxx mpi_cxx], [],
                 [has_mpicxx=no
                  MPICXX=
                  AC_MSG_WARN([
                  ----------------------------------------------------------
                   MPI C++ implementation does not support MPI-IO
                   Disabling C++ testing programs
                  ----------------------------------------------------------])]
          )])
          AC_LANG_POP(C++)
      fi
   fi
   if test "x${has_mpicxx}" = xyes ; then
      UD_CHECK_MPI_CPP_SEEK_SET
      UD_MSG_DEBUG(ac_cv_CHECK_MPI_CPP_SEEK_SET=$ac_cv_CHECK_MPI_CPP_SEEK_SET)
   fi
else
   dnl must invoke AC_PROG_CXX here otherwise got "error: conditional "am__fastdepCXX" was never defined."
   AC_PROG_CXX
   has_mpicxx=no
fi
AC_SUBST(has_mpicxx)
UD_MSG_DEBUG(has_mpicxx=$has_mpicxx)
AM_CONDITIONAL(HAS_MPICXX, [test x$has_mpicxx = xyes])

dnl I took the "-ansi" off the CFLAGS for the strict case, as it was
dnl preventing PATH_MAX from being defined for some of the test files,
dnl and I was having a tough time figuring out how to get it back! -- RobR
dnl Note this must be done after the type of C compiler is determined
AC_ARG_ENABLE(strict,
    [AS_HELP_STRING([--enable-strict],
                    [Turn on strict debugging with gcc. @<:@default: disabled@:>@])],
    [enable_strict=${enableval}], [enable_strict=no]
)

if test "x${enable_strict}" = xyes; then
   if test "x${GCC}" = xyes; then
      CFLAGS="${CFLAGS} -Wall -Wstrict-prototypes -Wmissing-prototypes -Wundef -Wpointer-arith -Wbad-function-cast"
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

if test "x${enable_fortran}" != xno ; then
   dnl if MPIF77 or F77 has not been set by users, then search from
   dnl CANDIDATE_MPIF77
   if test "x${MPIF77}" = x ; then
      UD_MPI_PATH_PROGS([MPIF77], [$CANDIDATE_MPIF77])
   fi
   dnl if MPIF90 or F90 has not been set by users, then search from
   dnl CANDIDATE_MPIF90
   if test "x${MPIF90}" = x ; then
      UD_MPI_PATH_PROGS([MPIF90], [$CANDIDATE_MPIF90])
   fi
fi
if test "x${enable_fortran}" = xyes ; then
   dnl user explicitly requests to build Fortran APIs
   if test "x${MPIF77}" = x && test "x${MPIF90}" = x ; then
      AC_MSG_ERROR([
         ------------------------------------------------------------
           Fortran support is explicitly requested, but configure
           cannot find a Fortran 77 or Fortran 90 compiler. Please
           specify the locations of the compilers either with the
           MPIF77 MPIF90 environment variables or the --with-mpi
           configure flag.
         ------------------------------------------------------------])
   fi
fi
UD_MSG_DEBUG(enable_fortran=${enable_fortran})

dnl Check if MPIF90 is a valid MPI compiler
mpi_mod=no
if test "x${MPIF90}" != x ; then
   FC=${MPIF90}
   AC_PROG_FC
   dnl FCFLAGS is set in AC_PROG_FC

   FCLIBS_save="$FCLIBS"
   FCLIBS=""
   AC_FC_LIBRARY_LDFLAGS
   UD_MSG_DEBUG([before FCLIBS=$FCLIBS])
   dnl The autoconf macro for finding FCLIBS sometimes makes mistakes
   dnl (particularly with the Fujitsu frt compiler).  This next step
   dnl first sees if the FCLIBS is valid with the Fortran compiler
   dnl This also happens to Solaris Studio Fortran compilers
   dnl AC_PROG_FC_FCLIBS_VALID
   dnl Now see if FCLIBS works with the C compiler (remove invalid ones)
   PAC_PROG_FC_CHECK_FCLIBS

   dnl replace FLIBS and F90LIBS with FCLIBS
   UD_MSG_DEBUG([after FCLIBS=$FCLIBS])
   FLIBS="$FCLIBS"
   dnl F90LIBS="$FCLIBS"

   AC_LANG_PUSH([Fortran])
   AC_MSG_CHECKING([if mpi.mod is available])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[        use mpi]])],
                     [mpi_mod=yes], [mpi_mod=no]
   )
   AC_MSG_RESULT($mpi_mod)
   AC_MSG_CHECKING([if MPIF90 is valid])
   if test "x${mpi_mod}" = xyes ; then
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[        use mpi
        integer err, rank
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)]])],
                           [valid_mpif90=yes],[valid_mpif90=no]
      )
   else
      dnl mpi.mod is not available, check if can use mpif.h
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[        include "mpif.h"
        integer err, rank
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)]])],
                           [valid_mpif90=yes],[valid_mpif90=no]
      )
   fi
   AC_MSG_RESULT($valid_mpif90)
   AC_LANG_POP([Fortran])
   if test "x${valid_mpif90}" = xno ; then
      if test "x${enable_fortran}" = xyes ; then
         AC_MSG_ERROR([
            ------------------------------------------------------------
              Invalid MPI Fortran 90 compiler specified: "${MPIF90}"
              A working MPI compiler is required. Please specify the
              location of one either with the MPIF90 environment
              variable or the --with-mpi configure flag.
            ------------------------------------------------------------])
      fi
      AC_MSG_WARN([
         ------------------------------------------------------------
           Invalid MPI Fortran 90 compiler: "${MPIF90}"
           Fortran feature is disabled
         ------------------------------------------------------------])
      enable_fortran=no
   else
      if test "x${mpi_mod}" = xyes ; then
         USE_MPIF_HEADER="use mpi, only: MPI_OFFSET_KIND"
      else
         USE_MPIF_HEADER="include \"mpif.h\""
      fi
      AC_SUBST(USE_MPIF_HEADER)
   fi
fi
AM_CONDITIONAL(HAVE_MPI_MOD, [test x$mpi_mod = xyes])

dnl Check if MPIF77 is a valid MPI compiler
ac_f77_support_freeform=no
if test "x${enable_fortran}" != xno && test "x${MPIF77}" != xno ; then
   F77=${MPIF77}
   AC_PROG_F77
   dnl FFLAGS is set in AC_PROG_F77

   AC_F77_LIBRARY_LDFLAGS

   FC_saved=${FC}
   FC=${F77}
   dnl customized AC_FC_FREEFORM: to just get ac_cv_fc_freeform without
   dnl appending ac_cv_fc_freeform to FCFLAGS
   dnl ac_cv_fc_freeform is the flag for enabling Fortran free form
   dnl we use this flag only in testing programs
   UD_FC_FREEFORM
   if test "x${ac_cv_fc_freeform}" != xunknown ; then
      ac_f77_support_freeform=yes
      FFREEFORMFLAG=${ac_cv_fc_freeform}
   fi
   UD_MSG_DEBUG([FC=$FC FFREEFORMFLAG=$FFREEFORMFLAG])
   AC_SUBST(FFREEFORMFLAG)
   FC=${FC_saved}

   AC_LANG_PUSH([Fortran 77])
   AC_MSG_CHECKING([if MPIF77 is valid])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[        include "mpif.h"
        integer comm
        comm = MPI_COMM_WORLD]])],
                     [valid_mpif77=yes], [valid_mpif77=no]
   )
   AC_MSG_RESULT($valid_mpif77)
   AC_LANG_POP([Fortran 77])
   if test "x${valid_mpif77}" = xno ; then
      if test "x${enable_fortran}" = xyes ; then
         AC_MSG_ERROR([
            ------------------------------------------------------------
              Invalid MPI Fortran 77 compiler specified: "${MPIF77}"
              A working MPI compiler is required. Please specify the
              location of one either with the MPIF77 environment
              variable or the --with-mpi configure flag
            ------------------------------------------------------------])
      fi
      AC_MSG_WARN([
         ------------------------------------------------------------
           MPI Fortran header file mpif.h is not valid:
           Fortran feature is disabled
         ------------------------------------------------------------])
      enable_fortran=no
   fi
fi
AM_CONDITIONAL(HAVE_F77_SUPPORT_FREEFORM, [test x$ac_f77_support_freeform = xyes])

if test "x${enable_fortran}" = xauto ; then
   if test "x${valid_mpif77}" = xyes && test "x${valid_mpif90}" = xyes ; then
      enable_fortran=yes
   else
      enable_fortran=no
      AC_MSG_WARN([
         ------------------------------------------------------------
           Either MPI Fortran 77 or 90 compiler is not working:
               "MPIF77 = ${MPIF77}"
               "MPIF90 = ${MPIF90}"
           Disable Fortran feature ...
         ------------------------------------------------------------])
   fi
fi

if test "x${enable_fortran}" = xyes ; then
   dnl GNU Fortran compiler automatically invokes preprocessor for files with
   dnl extension .F and .F90. To manually invoke the preprocessor on any file,
   dnl use compiler flag -cpp. To disable, use -nocpp.
   dnl For other compilers, we need to find the compile flags for Fortran
   dnl preprocessor.
   dnl AC_FC_PP_SRCEXT is first introduced in autoconf 2.69
   dnl steal AC_FC_PP_SRCEXT from autoconf 2.69 to make UD_FC_PP_SRCEXT
   UD_FC_PP_SRCEXT([F])    dnl sets ac_cv_fc_pp_srcext_F
   UD_FC_PP_SRCEXT([F90])  dnl sets ac_cv_fc_pp_srcext_F90

   dnl compiler command-line define preprocessor flag, result in ${FC_DEFINE}
   dnl Not all Fortran compilers use -D
   dnl AC_FC_PP_DEFINE is first introduced in autoconf 2.69
   UD_FC_PP_DEFINE

   dnl check compiler flags for file extensions in .f .F .f90 .F90
   AC_LANG_PUSH([Fortran])
   AC_FC_SRCEXT([f])
   AC_FC_SRCEXT([F])
   AC_FC_SRCEXT([f90])
   AC_FC_SRCEXT([F90])
   AC_LANG_POP([Fortran])

   dnl customized AC_FC_FIXEDFORM: to just get ac_cv_fc_fixedform without
   dnl appending ac_cv_fc_fixedform to FCFLAGS
   dnl ac_cv_fc_fixedform flag is for enabling Fortran fixed form.
   dnl we use this flag only in testing programs
   UD_FC_FIXEDFORM
   FFIXEDFORMFLAG=${ac_cv_fc_fixedform}
   UD_MSG_DEBUG([FC=$FC FFIXEDFORMFLAG=$FFIXEDFORMFLAG])
   AC_SUBST(FFIXEDFORMFLAG)

   # Checking for Fortran types also determines the Fortran name mangling
   # and places the value into FCALLSCSUB as the C name corresponding
   # to the Fortran name SUB 
   AC_FC_FUNC(sub, [FCALLSCSUB])

   # determine the correct name mapping
   case $FCALLSCSUB in
       SUB)
           AC_DEFINE(F77_NAME_UPPER,,[Define if Fortran names are uppercase])
           ;;
       sub_)
           # This is the hard case.  Gcc uses one _ unless the name includes
           # an underscore, in which case it gets two trailing underscores.
           # Use essentially the same configure code that the original configure
           # used to determine SUB
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
   dnl UD_CHECK_PGF77
   dnl UD_MSG_DEBUG([ac_cv_fc_compiler_pgf77=$ac_cv_fc_compiler_pgf77])

   dnl Check if the Fortran compiler is an NAG
   UD_CHECK_FC_NAG
   if test "x${ac_cv_fc_compiler_nag}" = xyes ; then
       FFLAGS="${FFLAGS}  -DNAGFortran -mismatch"
      FCFLAGS="${FCFLAGS} -DNAGFortran -mismatch"
   fi

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

   AC_SUBST(MPIF77)
   AC_SUBST(MPIF90)
fi

# compute canonical system types
AC_CANONICAL_BUILD
AC_CANONICAL_HOST

dnl libtool v2.4.6 was released in 2015-02-15
dnl Travis CI only has v2.4.2
LT_PREREQ([2.4.2])
dnl LT_INIT([dlopen disable-shared])
dnl LT_INIT([dlopen])
LT_INIT([disable-shared])

AC_PROG_INSTALL

have_yacc_lex=no
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
AM_CONDITIONAL(HAVE_YACC_LEX, [test x$have_yacc_lex = xyes])

dnl AC_PROG_AWK
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_SUBST(SET_MAKE)
AC_PROG_EGREP

dnl AC_PROG_SED and AC_PROG_GREP are only available on autoconf 2.60 and later
AC_PROG_SED
AC_PROG_GREP

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

dnl AM_PROG_AR is first recognized/traced in autoconf 2.68
dnl Use it only if we want support for unusual archivers such as Microsoft
dnl lib.
dnl AM_PROG_AR
dnl UD_PROG_AR()
dnl UD_PROG_NM()
dnl We could use the PAC check for ranlib (it also makes sure that ranlib
dnl works, which is not always true, particularly when GNU tools are
dnl installed on a system that does not have (or need) ranlib
dnl libtoolize: `AC_PROG_RANLIB' is rendered obsolete by `LT_INIT'
dnl AC_PROG_RANLIB

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
AC_CHECK_FUNCS([strerror access unlink])

AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
                    [Enable PnetCDF internal debug mode. This also enables safe mode.
                     @<:@default: disabled@:>@])],
    [debug=${enableval}], [debug=no]
)

PNETCDF_DEBUG=0
if test "x${debug}" = xyes; then
   dnl malloc memory allocation tracing relies on tdelete and tsearch
   AC_CHECK_HEADERS([search.h])
   AC_CHECK_FUNCS([tsearch tdelete])

   dnl check required functions for enabling malloc tracing
   if (test "x${ac_cv_func_tsearch}" = xyes) &&
      (test "x${ac_cv_func_tdelete}" = xyes) ; then
      AC_DEFINE(PNC_MALLOC_TRACE)
   fi
   PNETCDF_DEBUG=1
fi
AC_SUBST(PNETCDF_DEBUG)
AM_CONDITIONAL(PNETCDF_DEBUG, [test x"$PNETCDF_DEBUG" = x1])

dnl check the size of MPI_Offset. PnetCDF requires it be 8 bytes.
AC_CHECK_TYPE([MPI_Offset], [], [], [#include <mpi.h>])
if test "x${ac_cv_type_MPI_Offset}" = xyes; then
   AC_CHECK_SIZEOF([MPI_Offset], [], [#include <mpi.h>])
else
   AC_MSG_ERROR([Unable to find type MPI_Offset in mpi.h])
fi
if test "$ac_cv_sizeof_MPI_Offset" -lt "8"; then
   AC_MSG_ERROR([
         ------------------------------------------------------------
           Building of Parallel NetCDF cannot continue due to the
           size of MPI_Offset being less than 8 bytes. Please use
           an MPI library that supports large files.
         ------------------------------------------------------------])
fi

AC_CHECK_SIZEOF([MPI_Aint], [], [#include <mpi.h>])
AM_CONDITIONAL(SIZEOF_MPI_AINT_IS_4, [test x$ac_cv_sizeof_MPI_Aint = x4])

dnl the nonblocking routines build up lists of requests with MPI_Type_struct.
dnl If MPI_Offset not the same size as MPI_Aint, the arrays passed around will
dnl get mangled. 
if test "$ac_cv_sizeof_MPI_Offset" -ne "$ac_cv_sizeof_MPI_Aint"; then
   AC_MSG_WARN([
           MPI_Offset and MPI_Aint are detected of different sizes.
           The request aggregation feature in non-blocking APIs is
           thus disabled.])
   enable_aggregation=no
else
   AC_DEFINE(ENABLE_REQ_AGGREGATION,,[Define if able to support request aggregation in nonblocking routines])
   enable_aggregation=yes
fi

dnl check MPI-2 only functions
AC_CHECK_FUNCS([MPI_Info_dup \
                MPI_Info_free \
                MPI_Get_address \
                MPI_Type_create_subarray \
                MPI_Type_create_hvector \
                MPI_Type_create_hindexed \
                MPI_Type_create_struct \
                MPI_Type_create_resized \
                MPI_Type_get_extent])

dnl Check presence of MPI COMBINERS. These are of type int.
AC_CHECK_DECLS([MPI_COMBINER_DUP,
                MPI_COMBINER_HVECTOR_INTEGER,
                MPI_COMBINER_HINDEXED_INTEGER,
                MPI_COMBINER_SUBARRAY,
                MPI_COMBINER_DARRAY,
                MPI_COMBINER_RESIZED,
                MPI_COMBINER_STRUCT_INTEGER,
                MPI_COMBINER_INDEXED_BLOCK,
                MPI_COMBINER_F90_REAL,
                MPI_COMBINER_F90_INTEGER,
                MPI_COMBINER_F90_COMPLEX],
                [], [], [[#include <mpi.h>]])

dnl Check presence of various MPI error classes.
dnl These could be enums, so we have to do compile checks.
AC_CHECK_DECLS([MPI_ERR_FILE_EXISTS,
                MPI_ERR_NO_SUCH_FILE,
                MPI_ERR_AMODE,
                MPI_ERR_NOT_SAME,
                MPI_ERR_BAD_FILE,
                MPI_ERR_READ_ONLY,
                MPI_ERR_ACCESS,
                MPI_ERR_NO_SPACE,
                MPI_ERR_QUOTA],
                [], [], [[#include <mpi.h>]])

dnl Check presence of C MPI data types.
dnl Note they are of type MPI_Datatype, thus cannot use AC_CHECK_DECLS
UD_CHECK_MPI_DATATYPE(MPI_CHAR)
UD_CHECK_MPI_DATATYPE(MPI_BYTE)
UD_CHECK_MPI_DATATYPE(MPI_SIGNED_CHAR)
UD_CHECK_MPI_DATATYPE(MPI_UNSIGNED_CHAR)
UD_CHECK_MPI_DATATYPE(MPI_SHORT)
UD_CHECK_MPI_DATATYPE(MPI_UNSIGNED_SHORT)
UD_CHECK_MPI_DATATYPE(MPI_INT)
UD_CHECK_MPI_DATATYPE(MPI_UNSIGNED)
UD_CHECK_MPI_DATATYPE(MPI_LONG)
UD_CHECK_MPI_DATATYPE(MPI_FLOAT)
UD_CHECK_MPI_DATATYPE(MPI_DOUBLE)
UD_CHECK_MPI_DATATYPE(MPI_LONG_LONG_INT)
UD_CHECK_MPI_DATATYPE(MPI_UNSIGNED_LONG_LONG)
UD_CHECK_MPI_DATATYPE(MPI_UB)
UD_CHECK_MPI_DATATYPE(MPI_LB)
UD_CHECK_MPI_DATATYPE(MPI_OFFSET)
AM_CONDITIONAL(DECL_MPI_OFFSET, [test x$ac_cv_CHECK_MPI_DATATYPE_MPI_OFFSET = xyes])

dnl Check for presence of Fortran types
dnl These could be enums, so we have to do compile checks.
dnl
dnl We do this for a couple of reasons.  First, the MPI might have been
dnl built without Fortran support, in which case these types might not
dnl exist.  Second, we need to map these types to corresponding C types
dnl where possible to simplify processing at run time.
if test "x${enable_fortran}" = xyes ; then
    UD_CHECK_MPI_DATATYPE(MPI_CHARACTER)
    UD_CHECK_MPI_DATATYPE(MPI_REAL)
    UD_CHECK_MPI_DATATYPE(MPI_INTEGER)
    UD_CHECK_MPI_DATATYPE(MPI_DOUBLE_PRECISION)
    UD_CHECK_MPI_DATATYPE(MPI_INTEGER1)
    UD_CHECK_MPI_DATATYPE(MPI_INTEGER2)
    UD_CHECK_MPI_DATATYPE(MPI_INTEGER4)
    UD_CHECK_MPI_DATATYPE(MPI_INTEGER8)
    UD_CHECK_MPI_DATATYPE(MPI_INTEGER16)
    UD_CHECK_MPI_DATATYPE(MPI_REAL4)
    UD_CHECK_MPI_DATATYPE(MPI_REAL8)
    UD_CHECK_MPI_DATATYPE(MPI_REAL16)
    UD_CHECK_MPI_DATATYPE(MPI_COMPLEX8)
    UD_CHECK_MPI_DATATYPE(MPI_COMPLEX16)
    UD_CHECK_MPI_DATATYPE(MPI_COMPLEX32)
fi

AC_C_CHAR_UNSIGNED
AC_C_BIGENDIAN
AM_CONDITIONAL(IS_BIGENDIAN, [test x$ac_cv_c_bigendian = xyes])
AC_SUBST(ac_cv_c_bigendian)

in_place_swap=yes
AC_ARG_ENABLE([in-place-swap],
    [AS_HELP_STRING([--disable-in-place-swap],
                    [Disable memory in-place byte swap on Little Endian
                     machines. @<:@default: enabled@:>@])],
    [in_place_swap=${enableval}], [in_place_swap=yes]
)
if test "x${in_place_swap}" = xno ; then
   if test "x${ac_cv_c_bigendian}" = xyes ; then
      AC_MSG_WARN([--disable-in-place-swap takes no effect on Big Endian])
      in_place_swap=yes
   else
      AC_DEFINE(DISABLE_IN_PLACE_SWAP)
   fi
fi
AC_SUBST(in_place_swap)

dnl For big Endian, put buffer needs no byte swap and hence can be declared as
dnl INTENT(IN). For little Endian, put buffer may be used for byte swap in
dnl place and hence must be declared as INTENT(INOUT).
dnl This will configure/produce the file src/binding/f90/api.f90
if test "x${ac_cv_c_bigendian}" = xyes  || (test "x${in_place_swap}" = xno) ; then
   INTENTV="IN"
else
   INTENTV="INOUT"
fi
AC_SUBST(INTENTV)

AC_TYPE_SIZE_T
AC_TYPE_OFF_T
AC_CHECK_TYPES([ssize_t, ptrdiff_t, schar, uchar, ushort, uint, longlong, ulonglong, int64, uint64])

AC_CHECK_SIZEOF(size_t)
AC_CHECK_SIZEOF(off_t)
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

if test "x${enable_fortran}" = xyes ; then
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

    NFMPI_OFFSET="integer*$ac_cv_sizeof_MPI_Offset"
    AC_MSG_CHECKING([for Fortran NFMPI_OFFSET "$NFMPI_OFFSET"])
    dnl "
    AC_LANG_PUSH([Fortran 77])
    AC_COMPILE_IFELSE(
        [AC_LANG_SOURCE([
            subroutine sub(value)
            $NFMPI_OFFSET value
            end
        ])],
        [ac_cv_NFMPI_OFFSET=yes], [ac_cv_NFMPI_OFFSET=no]
    )
    AC_LANG_POP([Fortran 77])
    if test "$ac_cv_NFMPI_OFFSET" = yes ; then
        AC_MSG_RESULT(yes)
    else
        AC_MSG_RESULT(no)
        AC_MSG_ERROR([F77 does not support "$NFMPI_OFFSET"])
        dnl "
    fi
    ${RM} -rf conftest*
fi

SIZEOF_MPI_OFFSET=$ac_cv_sizeof_MPI_Offset
AC_SUBST(SIZEOF_MPI_OFFSET)

if test "x${enable_fortran}" = xyes ; then

   dnl AC_FC_MODULE_EXTENSION
   UD_FC_MODULE_EXTENSION
   dnl UD_FC_MODULE_EXTENSION defines FC_MODEXT
   if test "x${FC_MODEXT}" = x ; then
      AC_MSG_ERROR([cannot determine Fortran module file extension!])
   fi
   UD_MSG_DEBUG([FC_MODEXT=$FC_MODEXT])

   UD_FC_MODULE_FLAG
   dnl UD_FC_MODULE_FLAG defines FC_MODINC
   UD_MSG_DEBUG([FC_MODINC=$FC_MODINC])

   dnl AC_FC_MODULE_FLAG was first introduced in autoconf 2.69
   dnl AC_FC_MODULE_FLAG
   UD_FC_MODULE_OUTPUT_FLAG
   dnl UD_FC_MODULE_OUTPUT_FLAG defines FC_MODOUT
   dnl UD_MSG_DEBUG([FC_MODOUT=$FC_MODOUT])
fi

HAVE_F77_GNU_INT=no
HAVE_F77_INT1=no
HAVE_F77_INT2=no
HAVE_F77_INT8=no
if test "x${enable_fortran}" = xyes ; then
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

dnl
dnl Below is to check if a Fortran compiler produces module files with upper
dnl case file name, e.g. PNETCDF.mod. However, this does not work for Mac OSX
dnl file system which is case insensitive
dnl
ac_cv_prog_f90_uppercase_mod=no
if test "x${enable_fortran}" = xyes ; then
   UD_PROG_FC_UPPERCASE_MOD
fi
AM_CONDITIONAL(UPPER_CASE_MOD, [test x$ac_cv_prog_f90_uppercase_mod = xyes])

if test "x${enable_fortran}" = xyes ; then
   if test "x${enable_strict}" = xyes && test "x${GFC}" = xyes; then
      dnl not all Fortran compilers recognize -Wall
      FFLAGS="${FFLAGS} -Wall"
      FCFLAGS="${FCFLAGS} -Wall"
   fi
   AC_SUBST(FLIBS)
   AC_SUBST(FCLIBS)
   dnl AC_SUBST(F90LIBS)
   AC_SUBST(FLDFLAGS)
   dnl AC_SUBST(F90LDFLAGS)
fi

has_fortran=${enable_fortran}
AC_SUBST(has_fortran)
AM_CONDITIONAL(HAS_FORTRAN, [test x$has_fortran = xyes])

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
if test "x${enable_fortran}" = xyes ; then
   enable_f77=yes
   enable_fc=yes
fi
PAC_ENABLE_COVERAGE

if test "x${has_mpicxx}" = xyes ; then
   UD_CXX_MACRO_FUNC
   if test "x${ac_cv_cxx_macro_func}" == xyes ; then
      AC_DEFINE(HAVE_FUNC_MACRO)
   fi
   if test "x${ac_cv_cxx_macro_function}" == xyes ; then
      AC_DEFINE(HAVE_FUNCTION_MACRO)
   fi
fi

if test "x${debug}" = xyes; then
   dnl add -g flag if not presented
   dnl remove all -O and -fast flags
   dnl add -O0 to all flags
   if ! echo "${CFLAGS}" | ${EGREP} -q -- "-g" ; then
      CFLAGS="${CFLAGS} -g"
   fi
   CFLAGS=`echo $CFLAGS | sed 's/-O. *//g' | sed 's/-fast *//g'`
   CFLAGS="${CFLAGS} -O0"

   if test "x${has_mpicxx}" = xyes ; then
      if ! echo "${CXXFLAGS}" | ${EGREP} -q -- "-g" ; then
         CXXFLAGS="${CXXFLAGS} -g"
      fi
      CXXFLAGS=`echo $CXXFLAGS | sed 's/-O. *//g' | sed 's/-fast *//g'`
      CXXFLAGS="${CXXFLAGS} -O0"
   fi

   if test "x${enable_fortran}" = xyes ; then
      if ! echo "${FFLAGS}" | ${EGREP} -q -- "-g" ; then
         FFLAGS="${FFLAGS} -g"
      fi
      if ! echo "${FCFLAGS}" | ${EGREP} -q -- "-g" ; then
         FCFLAGS="${FCFLAGS} -g"
      fi

       FFLAGS=`echo $FFLAGS  | sed 's/-O. *//g' | sed 's/-fast *//g'`
      FCFLAGS=`echo $FCFLAGS | sed 's/-O. *//g' | sed 's/-fast *//g'`
       FFLAGS="${FFLAGS}  -O0"
      FCFLAGS="${FCFLAGS} -O0"
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

AC_ARG_ENABLE([erange-fill],
   [AS_HELP_STRING([--disable-erange-fill],
                   [Disable use of fill value when out-of-range type
                    conversion causes NC_ERANGE error. @<:@default: enabled@:>@])],
   [enable_erange_fill=${enableval}], [enable_erange_fill=yes]
)
ENABLE_ERANGE_FILL=0
if test "x$enable_erange_fill" = "xyes" ; then
   ENABLE_ERANGE_FILL=1
fi
AC_SUBST(ENABLE_ERANGE_FILL)
AM_CONDITIONAL(ENABLE_ERANGE_FILL, [test x$enable_erange_fill = xyes])

AC_ARG_ENABLE([relax-coord-bound],
   [AS_HELP_STRING([--enable-relax-coord-bound],
                   [Enable relaxed error NC_EINVALCOORDS to allow coordinate
                    start argument equal to dimension size when argument
                    count is zero. @<:@default: disabled@:>@])],
   [enable_relax_coord_bound=${enableval}], [enable_relax_coord_bound=no]
)
RELAX_COORD_BOUND=0
if test "x$enable_relax_coord_bound" = "xyes" ; then
   AC_DEFINE(RELAX_COORD_BOUND)
   RELAX_COORD_BOUND=1
fi
AC_SUBST(RELAX_COORD_BOUND)
AM_CONDITIONAL([RELAX_COORD_BOUND], [test "x$RELAX_COORD_BOUND" = x1])

AC_ARG_ENABLE([doxygen],
  [AS_HELP_STRING([--enable-doxygen],
    [Enable generation of documentation.])])
test "x$enable_doxygen" = xyes || enable_doxygen=no
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

AC_MSG_CHECKING([whether test should be built and run])
AC_ARG_ENABLE([testsets],
              [AS_HELP_STRING([--disable-testsets],
                              [do not build or run netCDF tests])])
test "x$enable_testsets" = xno || enable_testsets=yes
build_tests=$enable_testsets
AC_MSG_RESULT($build_tests)
AM_CONDITIONAL(BUILD_TESTSETS, [test x$build_tests = xyes])

AC_ARG_ENABLE([large-file-test],
    [AS_HELP_STRING([--enable-large-file-test],
                    [Enable testing for large (>4GB) file/variable I/O. Note
                     "make testing" can run very slow. @<:@default: disabled@:>@])],
    [large_file_test=${enableval}], [large_file_test=no]
)
AM_CONDITIONAL(BUILD_LARGE_FILE_TEST, [test x$large_file_test = xyes])

dnl build benchmark programs from within PnetCDF
AM_CONDITIONAL(BUILD_BENCHMARKS_IN_PNETCDF, [test yes = yes])

AC_ARG_VAR(TESTSEQRUN, [Run command (on one MPI process) for "make check" on cross-compile environment. Example: "aprun -n 1". @<:@default: none@:>@])
AC_SUBST(TESTSEQRUN)

AC_ARG_VAR(TESTMPIRUN, [MPI run command for "make ptest", @<:@default: mpiexec -n NP@:>@])
if test "x${TESTMPIRUN}" = x ; then
   dnl set default to "mpiexec -n NP"
   UD_MPI_PATH_PROGS([TESTMPIRUN], [mpiexec mpirun])
   TESTMPIRUN="${TESTMPIRUN} -n NP"
fi
TESTMPIRUN2=`echo $TESTMPIRUN | sed 's/ NP/ 2/g'`
TESTMPIRUN3=`echo $TESTMPIRUN | sed 's/ NP/ 3/g'`
TESTMPIRUN4=`echo $TESTMPIRUN | sed 's/ NP/ 4/g'`
TESTMPIRUN6=`echo $TESTMPIRUN | sed 's/ NP/ 6/g'`
TESTMPIRUN8=`echo $TESTMPIRUN | sed 's/ NP/ 8/g'`
TESTMPIRUN10=`echo $TESTMPIRUN | sed 's/ NP/ 10/g'`
AC_SUBST(TESTMPIRUN2)
AC_SUBST(TESTMPIRUN3)
AC_SUBST(TESTMPIRUN4)
AC_SUBST(TESTMPIRUN6)
AC_SUBST(TESTMPIRUN8)
AC_SUBST(TESTMPIRUN10)

AC_ARG_VAR(TESTOUTDIR, [Output file directory for "make check" and "make ptest", @<:@default: ./@:>@])
if test "x${TESTOUTDIR}" = x ; then
   dnl set default to current directory
   TESTOUTDIR=.
fi
AC_SUBST(TESTOUTDIR)

dnl find if gcc is available for compiling ncoffsets to run in sequential
AC_PATH_PROG([SEQ_CC], [gcc], [$MPICC])
AC_SUBST(SEQ_CC)

# Configuration Date
if test "x$SOURCE_DATE_EPOCH" != "x" ; then
    AC_SUBST([CONFIG_DATE]) CONFIG_DATE="`date -u -d "${SOURCE_DATE_EPOCH}"`"
else
    AC_SUBST([CONFIG_DATE]) CONFIG_DATE="`date`"
fi

dnl AC_CONFIG_HEADERS([src/binding/f77/nfconfig_inc])
AC_CONFIG_FILES(Makefile \
                libpnetcdf.pc \
                src/Makefile \
                src/dispatchers/Makefile \
                src/libs/Makefile \
                src/include/Makefile \
                src/include/pnetcdf.h \
                src/drivers/Makefile \
                src/drivers/ncmpi/Makefile \
                src/binding/Makefile \
                src/binding/cxx/Makefile \
                src/binding/f77/Makefile \
                src/binding/f77/pnetcdf.inc \
                src/binding/f90/Makefile \
                src/binding/f90/pnetcdf.f90 \
                src/binding/f90/api.f90 \
                src/binding/f90/nfmpi_constants.f90 \
                src/utils/Makefile \
                src/utils/ncoffsets/Makefile \
                src/utils/pnetcdf_version/Makefile \
                src/utils/ncmpivalid/Makefile \
                src/utils/ncmpidiff/Makefile \
                src/utils/ncmpidump/Makefile \
                src/utils/ncmpigen/Makefile \
                examples/Makefile \
                examples/C/Makefile \
                examples/tutorial/Makefile \
                examples/CXX/Makefile \
                examples/F77/Makefile \
                examples/F90/Makefile \
                benchmarks/Makefile \
                benchmarks/C/Makefile \
                benchmarks/FLASH-IO/Makefile \
                doc/Makefile \
                man/Makefile \
                test/Makefile \
                test/common/Makefile \
                test/C/Makefile \
                test/fandc/Makefile \
                test/nc_test/Makefile \
                test/cdf_format/Makefile \
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
)

# The following dependency is for configure.in and configure
# See autoconf manual 2.69, Section 4.8.5 Automatic Remaking
# AC_CONFIG_FILES([stamp-h], [echo timestamp > stamp-h])

dnl add executable permission to pnetcdf-config and check_install
AC_CONFIG_FILES([pnetcdf-config], [chmod u+x pnetcdf-config])
AC_CONFIG_FILES([check_install],  [chmod u+x check_install])

AC_OUTPUT
echo "------------------------------------------------------------------------------"

if test "x${enable_mpi_io_test}" = xno ; then
    AC_MSG_WARN([
       NOTE: disabling the MPI-IO test is a VERY bad idea.
       Please make sure you know what you are doing])
fi

echo \
"
   ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

   Features:  Build static libraries                  - ${enable_static}
              Build shared libraries                  - ${enable_shared}
              Build Fortran APIs                      - ${enable_fortran}
              Build C++ APIs                          - ${has_mpicxx}"
if test "x${enable_aggreg}" = xno; then
   echo "\
              Request aggregation in nonblocking APIs - OFF"
fi
if test "x${enable_erange_fill}" = xno; then
   echo "\
              Fill variables when NC_ERANGE occurs    - OFF"
fi
if test "x${enable_subfiling}" = xyes; then
   echo "\
              Build subfiling support                 - ON"
fi
if test "x${enable_relax_coord_bound}" = xyes; then
   echo "\
              Relax start coordinate bound check      - ON"
fi
if test "x${ac_cv_c_bigendian}" = xno  && (test "x${in_place_swap}" = xno) ; then
   echo "\
              Memory in-place byte swap               - OFF"
fi
if test "x${large_file_test}" = xyes; then
   echo "\
              Testing large file/variable I/O         - ON"
fi
if test "x${debug}" = xyes; then
   echo "\
              PnetCDF internal debug mode             - ON"
fi
if test "x${enable_fortran}" = xyes && (test "x${ac_f77_support_freeform}" = xno) ; then
   echo "\
              Support free form in Fortran 77         - OFF"
fi

echo "\

   Compilers: MPICC    = ${MPICC}"
if test "${has_mpicxx}" = yes ; then
   echo "\
              MPICXX   = ${MPICXX}"
fi
if test "${enable_fortran}" = yes ; then
   echo "\
              MPIF77   = ${MPIF77}
              MPIF90   = ${MPIF90}"
fi
echo "\
              CFLAGS   = ${CFLAGS}"
if test "x${CPPFLAGS}" != x ; then
   echo "\
              CPPFLAGS = ${CPPFLAGS}"
fi
if test "${has_mpicxx}" = yes ; then
   echo "\
              CXXFLAGS = ${CXXFLAGS}"
fi
if test "${enable_fortran}" = yes ; then
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
echo "\

   Now run 'make' to build the library and utility tools.
   Then run 'make @<:@<target>@:>@' for testing and installation, where the
   optional <target> can be:
              check    - build and run sequential test programs
              ptest    - build and run parallel test programs on 4 MPI processes
              tests    - build all test programs (build only, no run)
              install  - install PnetCDF library in ${prefix}
------------------------------------------------------------------------------"

