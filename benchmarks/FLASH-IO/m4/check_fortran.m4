dnl $Id$
dnl user-defined macros (UD_) for PnetCDF configure


dnl steal from autoconf 2.69
AC_DEFUN([UD_FC_PP_SRCEXT],
[AC_LANG_PUSH(Fortran)dnl
AC_CACHE_CHECK([for Fortran flag to compile preprocessed .$1 files],
                ac_cv_fc_pp_srcext_$1,
[ac_ext=$1
ac_fcflags_pp_srcext_save=$ac_fcflags_srcext
ac_fcflags_srcext=
ac_cv_fc_pp_srcext_$1=unknown
case $ac_ext in #(
  [[fF]]77) ac_try=f77-cpp-input;; #(
  *) ac_try=f95-cpp-input;;
esac
for ac_flag in none -ftpp -fpp -Tf "-fpp -Tf" -xpp=fpp -Mpreprocess "-e Z" \
               -cpp -xpp=cpp -qsuffix=cpp=$1 "-x $ac_try" +cpp -Cpp; do
  test "x$ac_flag" != xnone && ac_fcflags_srcext="$ac_flag"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#if 0
#include <ac_nonexistent.h>
      choke me
#endif]])],
    [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#if 1
#include <ac_nonexistent.h>
      choke me
#endif]])],
       [],
       [ac_cv_fc_pp_srcext_$1=$ac_flag; break])])
done
rm -f conftest.$ac_objext conftest.$1
ac_fcflags_srcext=$ac_fcflags_pp_srcext_save
])
if test "x$ac_cv_fc_pp_srcext_$1" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Fortran could not compile preprocessed .$1 files])])
else
  ac_fc_srcext=$1
  if test "x$ac_cv_fc_pp_srcext_$1" = xnone; then
    ac_fcflags_srcext=""
    FCFLAGS_[]$1[]=""
  else
    ac_fcflags_srcext=$ac_cv_fc_pp_srcext_$1
    FCFLAGS_[]$1[]=$ac_cv_fc_pp_srcext_$1
  fi
  AC_SUBST(FCFLAGS_[]$1)
  $2
fi
AC_LANG_POP(Fortran)dnl
])# UD_FC_PP_SRCEXT

dnl steal from autoconf 2.69
AC_DEFUN([UD_FC_PP_DEFINE],
[AC_LANG_PUSH([Fortran])dnl
ac_fc_pp_define_srcext_save=$ac_fc_srcext
ac_ext_saved=$ac_ext
UD_FC_PP_SRCEXT([F])
ac_ext=F
AC_CACHE_CHECK([how to define symbols for preprocessed Fortran],
  [ac_cv_fc_pp_define],
[ac_fc_pp_define_srcext_save=$ac_fc_srcext
ac_cv_fc_pp_define=unknown
ac_fc_pp_define_FCFLAGS_save=$FCFLAGS
for ac_flag in -D -WF,-D -Wp,-D -Wc,-D
do
  FCFLAGS="$ac_fc_pp_define_FCFLAGS_save ${ac_flag}FOOBAR ${ac_flag}ZORK=42"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#ifndef FOOBAR
      choke me
#endif
#if ZORK != 42
      choke me
#endif]])],
    [ac_cv_fc_pp_define=$ac_flag])
  test x"$ac_cv_fc_pp_define" != xunknown && break
done
FCFLAGS=$ac_fc_pp_define_FCFLAGS_save
])
ac_fc_srcext=$ac_fc_pp_define_srcext_save
if test "x$ac_cv_fc_pp_define" = xunknown; then
  FC_DEFINE=
  m4_default([$2],
             [AC_MSG_ERROR([Fortran does not allow to define preprocessor symbols], 77)])
else
  FC_DEFINE=$ac_cv_fc_pp_define
  $1
fi
ac_ext=$ac_ext_saved
AC_SUBST([FC_DEFINE])dnl
AC_LANG_POP([Fortran])dnl
])


dnl steal from autoconf 2.69
# UD_FC_MODULE_FLAG([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ---------------------------------------------------------------------
# Find a flag to include Fortran 90 modules from another directory.
# If successful, run ACTION-IF-SUCCESS (defaults to nothing), otherwise
# run ACTION-IF-FAILURE (defaults to failing with an error message).
# The module flag is cached in the ac_cv_fc_module_flag variable.
# It may contain significant trailing whitespace.
#
# Known flags:
# gfortran: -Idir, -I dir (-M dir, -Mdir (deprecated), -Jdir for writing)
# g95: -I dir (-fmod=dir for writing)
# SUN: -Mdir, -M dir (-moddir=dir for writing;
#                     -Idir for includes is also searched)
# HP: -Idir, -I dir (+moddir=dir for writing)
# IBM: -Idir (-qmoddir=dir for writing)
# Intel: -Idir -I dir (-mod dir for writing)
# Absoft: -pdir
# Lahey: -mod dir
# Cray: -module dir, -p dir (-J dir for writing)
#       -e m is needed to enable writing .mod files at all
# Compaq: -Idir
# NAGWare: -I dir
# PathScale: -I dir  (but -module dir is looked at first)
# Portland: -module dir (first -module also names dir for writing)
# Fujitsu: -Am -Idir (-Mdir for writing is searched first, then '.', then -I)
#                    (-Am indicates how module information is saved)
AC_DEFUN([UD_FC_MODULE_FLAG],[
AC_CACHE_CHECK([Fortran 90 module inclusion flag], [ac_cv_fc_module_flag],
[AC_LANG_PUSH([Fortran])
ac_cv_fc_module_flag=unknown
mkdir conftest.dir
cd conftest.dir
AC_COMPILE_IFELSE([[
      module conftest_module
      contains
      subroutine conftest_routine
      write(*,'(a)') 'gotcha!'
      end subroutine
      end module]],
  [cd ..
   ac_fc_module_flag_FCFLAGS_save=$FCFLAGS
   # Flag ordering is significant for gfortran and Sun.
   for ac_flag in -I '-I ' -M '-M ' -p '-mod ' '-module ' '-Am -I'; do
     # Add the flag twice to prevent matching an output flag.
     FCFLAGS="$ac_fc_module_flag_FCFLAGS_save ${ac_flag}conftest.dir ${ac_flag}conftest.dir"
     AC_COMPILE_IFELSE([[
      program main
      use conftest_module
      call conftest_routine
      end program]],
       [ac_cv_fc_module_flag="$ac_flag"])
     if test "$ac_cv_fc_module_flag" != unknown; then
       break
     fi
   done
   FCFLAGS=$ac_fc_module_flag_FCFLAGS_save
])
rm -rf conftest.dir
AC_LANG_POP([Fortran])
])
if test "$ac_cv_fc_module_flag" != unknown; then
  FC_MODINC=$ac_cv_fc_module_flag
  $1
else
  FC_MODINC=
  m4_default([$2],
    [AC_MSG_ERROR([unable to find compiler flag for module search path])])
fi
AC_SUBST([FC_MODINC])
# Ensure trailing whitespace is preserved in a Makefile.
AC_SUBST([ac_empty], [""])
AC_CONFIG_COMMANDS_PRE([case $FC_MODINC in #(
  *\ ) FC_MODINC=$FC_MODINC'${ac_empty}' ;;
esac])dnl
])

# _ACX_FC_MISMATCH([ACTION-IF-SUCCESS],
#                  [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Finds the compiler flag needed to allow routines to be called with different
# argument types. The result is either "unknown", or the actual compiler flag
# required to downgrade consistency checking of procedure argument lists, which
# may be an empty string.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The flag is cached in the acx_cv_[]_AC_LANG_ABBREV[]_mismatch_flag variable.
#
# Known flags:
# NAGWare: -mismatch
#
AC_DEFUN([_ACX_FC_MISMATCH],
  [_AC_FORTRAN_ASSERT()dnl
   m4_pushdef([acx_cache_var], [acx_cv_[]_AC_LANG_ABBREV[]_mismatch_flag])dnl
   AC_MSG_CHECKING([for _AC_LANG compiler flag needed to allow routines to dnl
be called with different argument types])
   AC_CACHE_VAL([acx_cache_var],
     [acx_cache_var=unknown
      acx_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
      AC_LANG_CONFTEST([AC_LANG_PROGRAM([], [[      implicit none
      integer a
      real b
      character c
      call foo1(a)
      call foo1(b)
      call foo1(c)]])])
      for acx_flag in '' -mismatch; do
        _AC_LANG_PREFIX[]FLAGS="${acx_save_[]_AC_LANG_PREFIX[]FLAGS} $acx_flag"
        AC_COMPILE_IFELSE([], [acx_cache_var=$acx_flag])
        test "x$acx_cache_var" != xunknown && break
      done
      rm -f conftest.$ac_ext
      _AC_LANG_PREFIX[]FLAGS=$acx_save_[]_AC_LANG_PREFIX[]FLAGS])
   AS_IF([test -n "$acx_cache_var"],
     [AC_MSG_RESULT([$acx_cache_var])],
     [AC_MSG_RESULT([none needed])])
   AS_VAR_IF([acx_cache_var], [unknown], [m4_default([$2],
     [AC_MSG_FAILURE([unable to detect _AC_LANG compiler flag needed to dnl
allow routines to be called with different argument types])])], [$1])
   m4_popdef([acx_cache_var])])

# ACX_FC_MISMATCH([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
AC_DEFUN([ACX_FC_MISMATCH],
  [AC_LANG_PUSH([Fortran])_ACX_FC_MISMATCH($@)AC_LANG_POP([Fortran])])

# ACX_F77_MISMATCH([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
AC_DEFUN([ACX_F77_MISMATCH],
  [AC_LANG_PUSH([Fortran 77])_ACX_FC_MISMATCH($@)AC_LANG_POP([Fortran 77])])
