dnl
dnl Check for a Standard C compiler.  Prefer a native one over the
dnl GNU one to reduce the chance that the environment variable LIBS
dnl will have to be set to reference the GNU C runtime library.
dnl
AC_DEFUN([UD_PROG_CC],
[
    # Because we must have a C compiler, we treat an unset CC
    # the same as an empty CC.
    case "${CC}" in
	'')
	    case `uname` in
		ULTRIX)
		    # The native ULTRIX C compiler isn't standard.
		    ccs='gcc cc'
		    ;;
		*)
		    # xlc is before c89 because AIX's sizeof(long long)
		    # differs between the two.
		    #
		    ccs='xlc c89 acc cc gcc'
		    ;;
	    esac
	    for cc in $ccs; do
		AC_CHECK_PROG(CC, $cc, $cc)
		case "$CC" in
		    '') ;;
		    *)  break
			;;
		esac
	    done
	    case "${CC}" in
		'')
		    AC_MSG_ERROR("Could not find C compiler")
		    ;;
	    esac
	    ;;
    esac
    #
    # On some systems, a discovered compiler nevertheless won't
    # work (due to licensing, for example); thus, we check the
    # compiler with a test program.
    #
    AC_MSG_CHECKING([C compiler "$CC"])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[]])],
                      [AC_MSG_RESULT(works)],
                      [AC_MSG_RESULT(failed to compile test program)])
    AC_SUBST(CC)
    case "$CC" in
	*gcc*)
	    GCC=yes		# Expected by autoconf(1) macros
	    ;;
    esac
    case `uname -sr` in
	'HP-UX A.09'*)
	    AC_DEFINE(_HPUX_SOURCE)
	    ;;
    esac
])

dnl
dnl like AC_TYPE_LONG_DOUBLE, except checks for 'long long'
dnl
AC_DEFUN([UD_C_LONG_LONG],
   [AC_MSG_CHECKING(for long long)
    AC_CACHE_VAL(ac_cv_c_long_long,
    [if test "$GCC" = yes; then
        ac_cv_c_long_long=yes
     else
        AC_RUN_IFELSE([AC_LANG_SOURCE([[int main() {
            long long foo = 0;
            return(sizeof(long long) < sizeof(long)); }]])],
           [ac_cv_c_long_long=yes],
           [ac_cv_c_long_long=no],
           [:])
     fi
    ])
    AC_MSG_RESULT($ac_cv_c_long_long)
    if test $ac_cv_c_long_long = yes; then
       AC_DEFINE(HAVE_LONG_LONG)
    fi
   ]
)

dnl Check for utility for generating makefile dependencies.
dnl Should only be used at the UPC.
dnl
AC_DEFUN([UD_PROG_CC_MAKEDEPEND],
[
    AC_MSG_CHECKING(how to make dependencies)
    case `uname -s` in
	IRIX*|OSF1)
	    CC_MAKEDEPEND='cc -M'
	    ;;
	SunOS)
	    case `uname -r` in
		4*)
		    CC_MAKEDEPEND='cc -M'
		    ;;
		5*|*)
		    CC_MAKEDEPEND='cc -xM'
		    ;;
	    esac
	    ;;
	ULTRIX)
	    case `uname -m` in
		RISC)
		    CC_MAKEDEPEND='cc -M'
		    ;;
		VAX)	# Can't handle prototypes in netcdf.h
		    ;;
	    esac
	    ;;
	AIX)	# Writes to .u files rather than standard out
	    ;;
	HP-UX)	# Writes escaped newlines to standard error
	    ;;
    esac
    case "${CC_MAKEDEPEND}" in
	'')
	    CC_MAKEDEPEND=false
	    ;;
    esac
    AC_MSG_RESULT($CC_MAKEDEPEND)
    AC_SUBST(CC_MAKEDEPEND)
])

# AX_C_FLOAT_WORDS_BIGENDIAN# added by:
#   Warren Turkal <wt@penguintechs.org>
#
# Copyright © 2006 Daniel Amelang <dan@amelang.net>
#
# Copying and distribution of this file, with or without modification, are
# permitted in any medium without royalty provided the copyright notice and
# this notice are preserved.
#
# This macro will detect if double variables are words packed in big endian
# order while the bits in the words are arranged in little endian order. This
# macro was added to support the ARM architecture. The FLOAT_WORDS_BIGENDIAN
# macro will be set to 1 if the word order is big endian. If the word order is
# not big endian, FLOAT_WORDS_BIGENDIAN will be not be set.
AC_DEFUN([AX_C_FLOAT_WORDS_BIGENDIAN],
  [AC_CACHE_CHECK(whether float word ordering is bigendian,
                  ax_cv_c_float_words_bigendian, [

ax_cv_c_float_words_bigendian=unknown
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[

double d = 9090423496703681033747047890550501147621169273561563201479712084405348
886581669527372346909785805625751702019124748742951693213050356065000232756451757
0778480236724525140520121371739201496540132640109977779420565776568942592.0;

]])], [

if grep noonsees conftest.$ac_objext >/dev/null ; then
  ax_cv_c_float_words_bigendian=yes
fi
if grep seesnoon conftest.$ac_objext >/dev/null ; then
  if test "$ax_cv_c_float_words_bigendian" = unknown; then
    ax_cv_c_float_words_bigendian=no
  else
    ax_cv_c_float_words_bigendian=unknown
  fi
fi

])])

case $ax_cv_c_float_words_bigendian in
  yes)
    m4_default([$1],
      [AC_DEFINE([FLOAT_WORDS_BIGENDIAN], 1,
                 [Define to 1 if your system stores words within floats
                  with the most significant word first])]) ;;
  no)
    $2 ;;
  *)
    m4_default([$3],
      [AC_MSG_ERROR([

Unknown float word ordering. You need to manually preset
ax_cv_c_float_words_bigendian=no (or yes) according to your system.

    ])]) ;;
esac

])# AX_C_FLOAT_WORDS_BIGENDIAN

dnl Check if C compiler is GNU
dnl According to gcc manual the command-line option to show version is --version
dnl
dnl % gcc --version
dnl gcc (Ubuntu 4.8.4-2ubuntu1~14.04.4) 4.8.4
dnl Copyright (C) 2013 Free Software Foundation, Inc.
dnl This is free software; see the source for copying conditions.  There is NO
dnl warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
dnl
AC_DEFUN([UD_CHECK_CC_IS_GCC],[
   if test "x$1" = x ; then
      TEST_CC=$CC
   else
      TEST_CC=$1
   fi
   AC_MSG_CHECKING([if compiler $TEST_CC is GNU gcc])
   ac_cv_cc_is_GCC=no
   ac_CC_VER=`$TEST_CC --version`
   ac_CC_VENDOR=`echo $ac_CC_VER | cut -s -d' ' -f1`
   if test "x${ac_CC_VENDOR}" = xgcc ; then
      ac_cv_cc_is_GCC=yes
   fi
   unset ac_CC_VER
   unset ac_CC_VENDOR
   unset TEST_CC
   AC_MSG_RESULT([$ac_cv_cc_is_GCC])
])

dnl Check if C compiler is CLANG
dnl According to clang manual the command-line option to show version is
dnl --version
dnl
dnl % clang --version
dnl Ubuntu clang version 3.4-1ubuntu3 (tags/RELEASE_34/final) (based on LLVM 3.4)
dnl Target: x86_64-pc-linux-gnu
dnl Thread model: posix
dnl
dnl or
dnl
dnl clang version 3.4.2 (tags/RELEASE_34/dot2-final)
dnl Target: x86_64-redhat-linux-gnu
dnl Thread model: posix
dnl
AC_DEFUN([UD_CHECK_CC_IS_CLANG],[
   if test "x$1" = x ; then
      TEST_CC=$CC
   else
      TEST_CC=$1
   fi
   AC_MSG_CHECKING([if compiler $TEST_CC is Clang])
   ac_cv_cc_is_CLANG=no
   ac_CC_VER=`$TEST_CC --version`
   ac_CC_VENDOR=`echo $ac_CC_VER | ${GREP} -w clang`
   if test "x${ac_CC_VENDOR}" != x ; then
      ac_cv_cc_is_CLANG=yes
   fi
   unset ac_CC_VER
   unset ac_CC_VENDOR
   unset TEST_CC
   AC_MSG_RESULT([$ac_cv_cc_is_CLANG])
])

dnl Check if C compiler is Intel icc
dnl According to icc manual the command-line option to show version is --version
dnl
dnl % icc --version
dnl icc (ICC) 17.0.0 20160721
dnl Copyright (C) 1985-2016 Intel Corporation.  All rights reserved.
dnl
AC_DEFUN([UD_CHECK_CC_IS_ICC],[
   if test "x$1" = x ; then
      TEST_CC=$CC
   else
      TEST_CC=$1
   fi
   AC_MSG_CHECKING([if compiler $TEST_CC is Intel icc])
   ac_cv_cc_is_ICC=no
   ac_CC_VER=`$TEST_CC --version`
   ac_CC_VENDOR=`echo $ac_CC_VER | cut -s -d' ' -f1`
   if test "x${ac_CC_VENDOR}" = xicc ; then
      ac_cv_cc_is_ICC=yes
   fi
   unset ac_CC_VER
   unset ac_CC_VENDOR
   unset TEST_CC
   AC_MSG_RESULT([$ac_cv_cc_is_ICC])
])

dnl Check if C compiler is Fujitsu fccpx based
dnl According to mpifccpx manual the command-line option to show version is
dnl -showme
dnl
dnl % mpifccpx --showme
dnl /opt/FJSVtclang/GM-1.2.0-24/bin/fccpx -Kident_mpi -mt ...
dnl
AC_DEFUN([UD_CHECK_CC_IS_FCCPX],[
   if test "x$1" = x ; then
      TEST_CC=$CC
   else
      TEST_CC=$1
   fi
   AC_MSG_CHECKING([if compiler $TEST_CC is Fujitsu fccpx])
   ac_cv_cc_is_FCCPX=no
   ac_CC_VENDOR=`$TEST_CC --showme 2> /dev/null | cut -d' ' -f1 | xargs -r basename`
   if test "x${ac_CC_VENDOR}" = xfccpx ; then
      ac_cv_cc_is_FCCPX=yes
   fi
   unset ac_CC_VENDOR
   unset TEST_CC
   AC_MSG_RESULT([$ac_cv_cc_is_FCCPX])
])

dnl Check if C compiler is IBM XL based
dnl According to xlc manual the command-line option to show version is
dnl
dnl % xlc -qversion
dnl IBM XL C/C++ for Blue Gene, V12.1
dnl Version: 12.01.0000.0011
dnl
AC_DEFUN([UD_CHECK_CC_IS_XLC],[
   if test "x$1" = x ; then
      TEST_CC=$CC
   else
      TEST_CC=$1
   fi
   AC_MSG_CHECKING([if compiler $TEST_CC is IBM XLC])
   ac_cv_cc_is_XLC=no
   ac_CC_VER=`$TEST_CC -qversion >& conftest.ver`
   ac_CC_VENDOR=`head -c 6 conftest.ver`
   if test "x${ac_CC_VENDOR}" = "xIBM XL" ; then
      ac_cv_cc_is_XLC=yes
   fi
   ${RM} -f conftest.ver
   unset ac_CC_VER
   unset ac_CC_VENDOR
   unset TEST_CC
   AC_MSG_RESULT([$ac_cv_cc_is_XLC])
])

dnl Check if MPI compiler is pgcc based
dnl According to pgcc manual the command-line option to show version is -V
dnl
dnl % pgcc -V
dnl
dnl pgcc 16.9-0 64-bit target on x86-64 Linux -tp p7
dnl The Portland Group - PGI Compilers and Tools
dnl Copyright (c) 2016, NVIDIA CORPORATION.  All rights reserved.
dnl
AC_DEFUN([UD_CHECK_CC_IS_PGCC],[
   if test "x$1" = x ; then
      TEST_CC=$CC
   else
      TEST_CC=$1
   fi
   AC_MSG_CHECKING([if compiler $TEST_CC is PGI pgcc])
   ac_cv_cc_is_PGCC=no
   ac_CC_VER=`$TEST_CC -V -c 2> /dev/null`
   ac_CC_VENDOR=`echo $ac_CC_VER | cut -s -d' ' -f1`
   if test "x${ac_CC_VENDOR}" = xpgcc ; then
      ac_cv_cc_is_PGCC=yes
   fi
   unset ac_CC_VER
   unset ac_CC_VENDOR
   unset TEST_CC
   AC_MSG_RESULT([$ac_cv_cc_is_PGCC])
])

dnl Check if C compiler is Oracle Solaris Studio
dnl According to cc manual the command-line option to show version is -V
dnl
dnl % cc -V
dnl cc: Sun C 5.13 Linux_i386 2014/10/20
dnl
AC_DEFUN([UD_CHECK_CC_IS_SOLARIS],[
   if test "x$1" = x ; then
      TEST_CC=$CC
   else
      TEST_CC=$1
   fi
   AC_MSG_CHECKING([if compiler $TEST_CC is Solaris cc])
   ac_cv_cc_is_SOLARIS=no
   ac_CC_VER="$($TEST_CC -V 2>&1)"
   ac_CC_VENDOR=`echo $ac_CC_VER | cut -s -d' ' -f2`
   if test "x${ac_CC_VENDOR}" = xSun ; then
      ac_cv_cc_is_SOLARIS=yes
   fi
   unset ac_CC_VER
   unset ac_CC_VENDOR
   unset TEST_CC
   AC_MSG_RESULT([$ac_cv_cc_is_SOLARIS])
])

dnl Check C compiler base
dnl
AC_DEFUN([UD_CHECK_CC_BASE_VENDOR],[
   if test "x$1" = x ; then
      TEST_CC=$CC
   else
      TEST_CC=$1
   fi
   AC_MSG_CHECKING([compiler $TEST_CC base])
   ac_cv_cc_base_vendor=
   # Check GCC
   ac_CC_VER="$($TEST_CC --version 2>&1)"
   ac_CC_VER=`echo $ac_CC_VER | ${GREP} -w gcc`
   # AC_MSG_NOTICE(GCC ac_CC_VER=$ac_CC_VER)
   if test "x${ac_CC_VER}" != x ; then
      ac_cv_cc_base_vendor="GCC"
   else
      # Check CLANG
      ac_CC_VER="$($TEST_CC --version 2>&1)"
      ac_CC_VER=`echo $ac_CC_VER | ${GREP} -w clang`
      # AC_MSG_NOTICE(clang ac_CC_VER=$ac_CC_VER)
      if test "x${ac_CC_VER}" != x ; then
         ac_cv_cc_base_vendor="CLANG"
      else
         # Check Intel C
         ac_CC_VER="$($TEST_CC --version 2>&1)"
         # grep keyword Intel instead of "icc"
         ac_CC_VER=`echo $ac_CC_VER | ${GREP} -w Intel`
         # AC_MSG_NOTICE(icc ac_CC_VER=$ac_CC_VER)
         if test "x${ac_CC_VER}" != x ; then
            ac_cv_cc_base_vendor="ICC"
         else
            # Check XLC
            ac_CC_VER="$($TEST_CC -qversion 2>&1)"
            ac_CC_VER=`echo $ac_CC_VER | ${GREP} "IBM XL C"`
            # AC_MSG_NOTICE(XLC ac_CC_VER=$ac_CC_VER)
            if test "x${ac_CC_VER}" != x ; then
               ac_cv_cc_base_vendor="XLC"
            else
               # Check PGCC
               ac_CC_VER="$($TEST_CC -V -c 2>&1)"
               # grep keyword PGI instead of "pgcc"
               ac_CC_VER=`echo $ac_CC_VER | ${GREP} -w PGI`
               # AC_MSG_NOTICE(pgcc ac_CC_VER=$ac_CC_VER)
               if test "x${ac_CC_VER}" != x ; then
                  ac_cv_cc_base_vendor="PGCC"
               else
                  # Check SOLARIS
                  ac_CC_VER="$($TEST_CC -V 2>&1)"
                  ac_CC_VER=`echo $ac_CC_VER | ${GREP} -w Sun`
                  # AC_MSG_NOTICE(Sun ac_CC_VER=$ac_CC_VER)
                  if test "x${ac_CC_VER}" != x ; then
                     ac_cv_cc_base_vendor="SOLARIS"
                  else
                     # Check FCCPX
                     ac_CC_VER="$($TEST_CC --showme 2>&1)"
                     ac_CC_VER=`echo $ac_CC_VER | ${GREP} -w fccpx`
                     # AC_MSG_NOTICE(fccpx ac_CC_VER=$ac_CC_VER)
                     if test "x${ac_CC_VER}" != x ; then
                        ac_cv_cc_base_vendor="FCCPX"
                     else
                        # If just cc, check if it is a wrapper of GCC
                        ac_CC_VER="$($TEST_CC -v 2>&1)"
                        UD_MSG_DEBUG(GCC ac_CC_VER=$ac_CC_VER)
                        ac_CC_VER=`echo $ac_CC_VER | ${GREP} -w gcc`
                        if test "x${ac_CC_VER}" != x ; then
                           ac_cv_cc_base_vendor="GCC"
                        fi
                     fi
                  fi
               fi
            fi
         fi
      fi
   fi
   unset TEST_CC
   unset ac_CC_VER
   if test "x$ac_cv_cc_base_vendor" = x ; then
      AC_MSG_RESULT([unknown])
   else
      AC_MSG_RESULT([$ac_cv_cc_base_vendor])
   fi
])
