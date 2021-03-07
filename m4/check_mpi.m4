dnl
dnl check the availability of one MPI executable in $2
dnl
dnl $2 can be a single command, This is the case when user set the environment.
dnl The variable may contain the executable name followed by zeor or more
dnl command-line options. In the latter case, we check the first string token,
dnl the command name, and ignore the rest command-line options. For example,
dnl UD_MPI_PATH_PROG([MPICC], [mpicc -O2])
dnl
dnl In addition, the first token of $2 may contain the full path of
dnl the command. For example, UD_MPI_PATH_PROG([MPICC], [/usr/bin/mpicc -O2])
dnl
AC_DEFUN([UD_MPI_PATH_PROG], [
   if test "x$2" = x ; then
      AC_MSG_ERROR("2nd argument cannot be NULL")
   fi

   dnl 1st token in $2 must be the program name, rests are command-line options
   ac_first_token=`echo $2 | cut -d" " -f1`
   ac_rest_tokens=`echo $2 | cut -d" " -s -f2-`
   UD_MSG_DEBUG(ac_first_token=$ac_first_token) dnl executable name
   UD_MSG_DEBUG(ac_rest_tokens=$ac_rest_tokens) dnl command-line option

   dnl First check if ac_first_token contain a full path
   dnl If yes, check, check if the file exists. Need not check MPI_INSTALL.
   ac_mpi_prog_path=`AS_DIRNAME(["$ac_first_token"])`
   if test "x$ac_mpi_prog_path" != "x." ; then
      AC_MSG_CHECKING([whether $ac_first_token exists and is executable])
      if test -x "$ac_first_token" ; then
         AC_MSG_RESULT([yes])
         $1="$2"
      else
         AC_MSG_RESULT([no])
         $1=
      fi
   else
      dnl ac_first_token does not contain a full path
      ac_mpi_prog_$1=
      if test "x$MPI_INSTALL" != x ; then
         dnl First, check if it can be found under $MPI_INSTALL, i.e.
         dnl --with-mpi is used on configure command line
         if test -d "${MPI_INSTALL}/bin" ; then
            AC_MSG_CHECKING([$ac_first_token under ${MPI_INSTALL}/bin])
            if test -x "$MPI_INSTALL/bin/$ac_first_token" ; then
               AC_MSG_RESULT([yes])
               ac_mpi_prog_$1=$MPI_INSTALL/bin/$ac_first_token
            else
               AC_MSG_RESULT([no])
            fi
         else
            dnl ${MPI_INSTALL}/bin does not exist, search $MPI_INSTALL
            AC_MSG_CHECKING([$ac_first_token under ${MPI_INSTALL}])
            if test -x "$MPI_INSTALL/$ac_first_token" ; then
               AC_MSG_RESULT([yes])
               ac_mpi_prog_$1=$MPI_INSTALL/$ac_first_token
            else
               AC_MSG_RESULT([no])
            fi
         fi
         if test "x$ac_mpi_prog_$1" != x ; then
            if test "x$ac_rest_tokens" != x ; then
               $1="${ac_mpi_prog_$1} $ac_rest_tokens"
            else
               $1=${ac_mpi_prog_$1}
            fi
         else
            $1=
         fi
      else
         dnl MPI_INSTALL is not set, i.e. --with-mpi is not used
         AC_PATH_PROG([ac_mpi_prog_$1], [$ac_first_token])
         if test "x$ac_mpi_prog_$1" != x ; then
            if test "x$ac_rest_tokens" != x ; then
               $1="${ac_mpi_prog_$1} $ac_rest_tokens"
            else
               $1=${ac_mpi_prog_$1}
            fi
         else
            $1=
         fi
      fi
   fi
])

dnl
dnl check the availability of a list of MPI executables
dnl
dnl Note $2 can be a list of executable commands to be searched with each
dnl command being the executable file name without command-line option. This is
dnl the case when user does not set the environment variable, for example
dnl MPICC, and we must search one from the candidate list. For example,
dnl UD_MPI_PATH_PROGS([MPICC], [mpicc mpixlc mpifccpx mpipgcc])
dnl
AC_DEFUN([UD_MPI_PATH_PROGS], [
   ac_mpi_prog_$1=
   if test "x$MPI_INSTALL" != x ; then
      UD_MSG_DEBUG(--with-mpi=$MPI_INSTALL is used)
      if test -d "${MPI_INSTALL}/bin" ; then
         UD_MSG_DEBUG(search $2 under $MPI_INSTALL/bin)
         AC_PATH_PROGS([ac_mpi_prog_$1], [$2], [], [$MPI_INSTALL/bin])
      else
         dnl ${MPI_INSTALL}/bin does not exist, search $MPI_INSTALL
         UD_MSG_DEBUG(search $2 under $MPI_INSTALL)
         AC_PATH_PROGS([ac_mpi_prog_$1], [$2], [], [$MPI_INSTALL])
      fi
   else
      UD_MSG_DEBUG(--with-mpi=$MPI_INSTALL is NOT used)
      UD_MSG_DEBUG(search $2 under $PATH)
      AC_PATH_PROGS([ac_mpi_prog_$1], [$2])
   fi
   UD_MSG_DEBUG([ac_mpi_prog_$1=${ac_mpi_prog_$1}])
   if test "x${ac_mpi_prog_$1}" = x ; then
      dnl AC_CHECK_FILES fails when $2 is not found in cross compile
      dnl AC_CHECK_FILES([$2], [ac_mpi_prog_$1=$2])
      AC_PATH_PROGS([ac_mpi_prog_$1], [$2])
      dnl AC_CHECK_PROGS([ac_mpi_prog_$1], [$2])
      dnl AC_CHECK_PROGS([ac_mpi_prog_$1], [$2], [], [/])
      dnl ac_first_token=`echo $2 | cut -d" " -f1`
      dnl UD_MSG_DEBUG(check first token $ac_first_token of $2)
      dnl if test -f $ac_first_token ; then
         dnl UD_MSG_DEBUG(use file $ac_first_token as it exits)
         dnl ac_mpi_prog_$1=$2
      dnl fi
   fi
   $1=${ac_mpi_prog_$1}
])

dnl Check for presence of an MPI constant.
dnl These could be enums, so we have to do compile checks.
AC_DEFUN([UD_HAS_MPI_CONST], [
   AC_MSG_CHECKING(whether MPI constant $1 is defined )
   AC_COMPILE_IFELSE(
      [AC_LANG_SOURCE([
          #include <mpi.h>
          int dummy = $1;
      ])],
      [AC_MSG_RESULT(yes)
       AC_DEFINE(HAVE_$1, 1, available)
      ],
      [AC_MSG_RESULT(no)]
   )]
)

dnl Check for presence of an MPI datatype.
dnl These could be enums, so we have to do compile checks.
AC_DEFUN([UD_CHECK_MPI_DATATYPE], [
   AC_MSG_CHECKING(whether MPI datatype $1 is defined )
   AC_COMPILE_IFELSE(
      [AC_LANG_SOURCE([
          #include <mpi.h>
          MPI_Datatype dummy = $1;
      ])],
      [ac_cv_CHECK_MPI_DATATYPE_$1=yes],
      [ac_cv_CHECK_MPI_DATATYPE_$1=no]
   )
   AC_MSG_RESULT($ac_cv_CHECK_MPI_DATATYPE_$1)
   if test "x$ac_cv_CHECK_MPI_DATATYPE_$1" = xyes; then
       AC_DEFINE(HAVE_DECL_$1, 1, available)
   fi
])

dnl Check if older Intel MPI C compiler (4.x) for issue of redefined SEEK_SET
dnl See https://software.intel.com/en-us/articles/intel-cluster-toolkit-for-linux-error-when-compiling-c-aps-using-intel-mpi-library-compilation-driver-mpiicpc
AC_DEFUN([UD_CHECK_MPI_CPP_SEEK_SET], [
   AC_MSG_CHECKING(whether MPI C++ compiler redefines SEEK_SET )
   CXX=${MPICXX}
   AC_LANG_PUSH(C++)
   AC_COMPILE_IFELSE(
      [AC_LANG_SOURCE([
          #include <stdio.h>
          #include <mpi.h>
          int main() { return 0; }
      ])],
      [ac_cv_CHECK_MPI_CPP_SEEK_SET=no],
      [ac_cv_CHECK_MPI_CPP_SEEK_SET=yes]
   )
   AC_MSG_RESULT([$ac_cv_CHECK_MPI_CPP_SEEK_SET])
   AC_LANG_POP(C++)
])

#
# To check whether MPI C compiler supports shared libraries, below three
# functions LT_AC_CHECK_SHLIB, LT_AH_CHECK_SHLIB, and LT_AC_LINK_SHLIB_IFELSE
# are from http://lists.gnu.org/archive/html/libtool/2004-10/msg00222.html
#
# Usage example:
#
#  LT_AC_CHECK_SHLIB(mpi, MPI_Init, [], [AC_MSG_ERROR([
#    -----------------------------------------------------------------------
#      Building shared libraries is requested, but the MPI C compiler:
#      "${MPICC}"
#      is not built with support of shared libraries. Abort.
#    -----------------------------------------------------------------------])])
#

# LT_MPI_CHECK_SHLIB
# -----------------------------------------------------------------
# Try to link an MPI program using libtool. This function is useful for
# detecting whether the MPI library is built with shared library support.
AC_DEFUN([LT_MPI_CHECK_SHLIB],[
   AC_MSG_CHECKING([whether MPI library is built with shared library support])
   AC_LANG_CONFTEST([AC_LANG_PROGRAM([[#include <mpi.h>]], [[MPI_Init(0, 0);]])])
   $RM -rf $objdir conftest.$ac_objext conftest.la
   dnl RM must have -f option when calling libtool
   ac_RM_saved=${RM}
   if test "x$RM" = xrm || test "x$RM" = "x/bin/rm" ; then
      RM="$RM -f"
   fi
   ac_ltcompile='./libtool --mode=compile $MPICC -c $CFLAGS $CPPFLAGS conftest.$ac_ext -o conftest.lo >&AS_MESSAGE_LOG_FD'
   ac_ltlink_la='./libtool --mode=link $MPICC -rpath `pwd` $CFLAGS $LDFLAGS -o libconftest.la conftest.lo $LIBS >&AS_MESSAGE_LOG_FD'
   AS_IF([AC_TRY_EVAL([ac_ltcompile]) &&
       AC_TRY_EVAL([ac_ltlink_la]) &&
       AC_TRY_COMMAND([test -s libconftest.la])],
      [ac_cv_lt_mpi_check_shlib=yes],
      [ac_cv_lt_mpi_check_shlib=no])
   RM=$ac_RM_saved
   unset ac_RM_saved
   $RM -rf $objdir conftest* libconftest*
   AC_MSG_RESULT([$ac_cv_lt_mpi_check_shlib])
])# LT_MPI_CHECK_SHLIB

# CHECK_MPI_VERSION
# -----------------------------------------------------------------
# check MPI version and vendor
AC_DEFUN([CHECK_MPI_VERSION],[
   AC_REQUIRE([AX_COMPILER_VENDOR])
   AC_REQUIRE([AC_PROG_GREP])
   AC_MSG_CHECKING([MPI Standard version implemented])
   AC_COMPUTE_INT([mpi_version], [MPI_VERSION], [[#include <mpi.h>]])
   AC_COMPUTE_INT([mpi_subversion], [MPI_SUBVERSION], [[#include <mpi.h>]])
   if test "x$mpi_version" = x ; then
      AC_MSG_RESULT([information unavailable])
   else
      AC_MSG_RESULT([${mpi_version}.${mpi_subversion}])
   fi

   AC_CHECK_DECL([MPICH_VERSION],      [], [], [#include <mpi.h>])
   AC_CHECK_DECL([MPICH2_VERSION],     [], [], [#include <mpi.h>])
   AC_CHECK_DECL([OMPI_MAJOR_VERSION], [], [], [#include <mpi.h>])
   AC_CHECK_DECL([MVAPICH2_VERSION],   [], [], [#include <mpi.h>])
   AC_MSG_CHECKING([MPI vendor])

   # MACRO_FLAG is the CPP flag to show macro definitions. For GCC, Intel, and
   # PGI, it is -dM and prints to stdout. For Oracle Solaris Studio compiler,
   # it is -xdumpmacros=defs and prints to stderr
   MACRO_FLAG="-dM"
   if test "x$ax_cv_c_compiler_vendor" = xibm ; then
      MACRO_FLAG="-qshowmacros"
   fi

   saved_CPPFLAGS=$CPPFLAGS
   CPPFLAGS="$CPPFLAGS $MACRO_FLAG"
   AC_PREPROC_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]], [])],
                     [`cp conftest.i saved_conftest.i`])
   CPPFLAGS=$saved_CPPFLAGS
   unset MACRO_FLAG

   # version_str=`${EGREP} -E 'MVAPICH2_VERSION|MPICH_VERSION|MPICH2_VERSION|OMPI_MAJOR_VERSION|OMPI_MINOR_VERSION|OMPI_RELEASE_VERSION' conftest.i`
   # echo a variable deletes any linefeeds in that variable, so we cannot use
   # `echo $version_str | ${GREP} MVAPICH2_VERSION`. Instead, we can use
   # version=`${GREP} MPICH_VERSION <<< "$version_str" | cut -d' ' -d'"' -f2`

   # Note MVAPICH2's mpi.h also defines MPICH_VERSION, so this check must be
   # done before MPICH.
   if test -f saved_conftest.i ; then
      if test "x$ac_cv_have_decl_MVAPICH2_VERSION" = xyes ; then
         mvapich2_version=`${GREP} MVAPICH2_VERSION saved_conftest.i | cut -d' ' -d'"' -f2`
         AC_MSG_RESULT(MVAPICH2 $mvapich2_version)
         unset mvapich2_version
      elif test "x$ac_cv_have_decl_MPICH_VERSION" = xyes ; then
         mpich_version=`${GREP} MPICH_VERSION saved_conftest.i | cut -d' ' -d'"' -f2`
         AC_MSG_RESULT(MPICH $mpich_version)
         unset mpich_version
      elif test "x$ac_cv_have_decl_MPICH2_VERSION" = xyes ; then
         mpich2_version=`${GREP} MPICH2_VERSION saved_conftest.i | cut -d' ' -d'"' -f2`
         AC_MSG_RESULT(MPICH2 $mpich2_version)
         unset mpich2_version
      elif test "x$ac_cv_have_decl_OMPI_MAJOR_VERSION" = xyes ; then
         # AC_COMPUTE_INT([OMPI_MAJOR], [OMPI_MAJOR_VERSION], [[#include <mpi.h>]])
         # AC_COMPUTE_INT([OMPI_MINOR], [OMPI_MINOR_VERSION], [[#include <mpi.h>]])
         # AC_COMPUTE_INT([OMPI_RELEASE], [OMPI_RELEASE_VERSION], [[#include <mpi.h>]])
         OMPI_MAJOR=`${GREP} OMPI_MAJOR_VERSION saved_conftest.i | cut -d' ' -f3`
         OMPI_MINOR=`${GREP} OMPI_MINOR_VERSION saved_conftest.i | cut -d' ' -f3`
         OMPI_RELEASE=`${GREP} OMPI_RELEASE_VERSION saved_conftest.i | cut -d' ' -f3`
         ompi_version="${OMPI_MAJOR}.${OMPI_MINOR}.${OMPI_RELEASE}"
         AC_MSG_RESULT(OpenMPI $ompi_version)
         unset OMPI_MAJOR
         unset OMPI_MINOR
         unset OMPI_RELEASE
         unset ompi_version
      fi
   else
      AC_MSG_RESULT([unknown])
   fi
   ${RM} -f saved_conftest.i
])

dnl It is not sufficient to use AC_CHECK_DECLS to check whether deprecated MPI
dnl constants are defined. OpenMPI still define the constants deprecated by MPI
dnl 3.0 as some error messages. We need to call AC_COMPILE_IFELSE to check
dnl whether they can be used without compilation errors.
dnl
m4_define([_UD_CHECK_MPI_CONSTANTS],
   [AC_MSG_CHECKING([whether $1 is defined])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([$4],[[int dummy=$1;]])],
                      [ac_have_const=yes], [ac_have_const=no])
    AC_MSG_RESULT([$ac_have_const])
    if test "x$ac_have_const" = xyes ; then
       AC_DEFINE([HAVE_$1], [1], [Define if $1 is defined and not deprecated])
    fi]
   [m4_ifvaln([$2$3], [AS_IF([test x$ac_have_const = xyes], [$2], [$3])])])

AC_DEFUN([UD_CHECK_MPI_CONSTANTS],
   [m4_map_args_sep([_$0(], [, [$2], [$3], [$4])], [], $1)])

# MPI_COMPILER_BASE
# -----------------------------------------------------------------
# Returns the base compiler command used in MPI compiler wrapper in
# ac_cv_mpi_compiler_base_MPICC, ac_cv_mpi_compiler_base_MPICXX,
# ac_cv_mpi_compiler_base_MPIF90, such as /usr/bin/gcc, /usr/bin/g++,
# /usr/bin/gfortran
#
AC_DEFUN([MPI_COMPILER_BASE],[
   # before checking, remove compile command-line options, if there is any
   compile_cmd=`echo $$1 | cut -d" " -f1`
   AC_MSG_CHECKING([base compiler command in $1 wrapper])
   compile_basename=
   case $compile_cmd in
        mpicc | mpicxx | mpif77 | mpif90 | mpifort | *[[\\/]]mpicc | *[[\\/]]mpicxx | *[[\\/]]mpif77 | *[[\\/]]mpif90 | *[[\\/]]mpifort )
           # MPICH, OpenMPI
           compile_basename=`$compile_cmd -show | cut -d' ' -f1`
           ;;
        mpixlc | mpixlcxx | mpixlf77 | mpixlf90 | *[[\\/]]mpixlc | *[[\\/]]mpixlcxx | *[[\\/]]mpixlf77 | *[[\\/]]mpixlf90 )
           # IBM XL MPI compilers
           compile_basename=`$compile_cmd -show | cut -d' ' -f1`
           ;;
        mpifccpx | mpiFCCpx | mpifrtpx | *[[\\/]]mpifccpx | *[[\\/]]mpiFCCpx | *[[\\/]]mpifrtpx )
           # Fujitsu MPI compilers: fccpx, FCCpx, frtpx
           compile_basename=`$compile_cmd -showme | cut -d' ' -f1`
           ;;
        cc | CC | ftn | *[[\\/]]cc | *[[\\/]]CC | *[[\\/]]ftn )
           # For Cray PrgEnv-intel, cc is a wrapper of icc
           # For Cray PrgEnv-gnu, cc is a wrapper of gcc
           eval "$compile_cmd --version" < /dev/null >& conftest.ver
           compile_basename=`head -n1 conftest.ver |cut -d' ' -f1`
           ${RM} -f conftest.ver
           if test "x${compile_basename}" = x ; then
              # For Cray PrgEnv-cray, cc is a wrapper of Cray CC
              # Cray cc -V sends the output to stderr.
              eval "$compile_cmd -V" < /dev/null >& conftest.ver
              compile_basename=`head -n1 conftest.ver |cut -d' ' -f1`
              ${RM} -f conftest.ver
           fi
           ;;
        *) break;;
   esac
   if test "x${compile_basename}" != x ; then
      AC_MSG_RESULT([$compile_basename])
      AC_PATH_PROG([ac_cv_mpi_compiler_base_$1], [$compile_basename])
   else
      AC_MSG_RESULT([not found])
   fi
   unset compile_basename
   unset compile_cmd
])# MPI_COMPILER_BASE

