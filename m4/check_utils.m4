dnl $Id$
dnl UD macros for PnetCDF configure


dnl Convert a string to all uppercase.
dnl
define([uppercase],
[translit($1, abcdefghijklmnopqrstuvwxyz, ABCDEFGHIJKLMNOPQRSTUVWXYZ)])

dnl
dnl Check for an m4(1) preprocessor utility.
dnl
AC_DEFUN([UD_PROG_M4],[
   dnl AS_MESSAGE([checking for m4 preprocessor...])
   case "${M4-unset}" in
      unset) AC_CHECK_PROGS(M4, m4 gm4, []) ;;
      *) AC_CHECK_PROGS(M4, $M4 m4 gm4, []) ;;
   esac
   if test -z "$M4" ; then
      AC_MSG_ERROR("m4 utility program is required by PnetCDF")
   fi
   AC_MSG_CHECKING(m4 additional flags)
   case "${M4FLAGS-unset}" in
       unset) dnl test if M4 runs fine without option -B10000
              `${M4} /dev/null > conftest.err 2>&1`
              ac_cv_m4_stdout=`cat conftest.err`
              if test "x$ac_cv_m4_stdout" != x; then
                 M4FLAGS=-B10000
              fi
              ${RM} -f conftest.err
              ;;
   esac
   if test "x$M4FLAGS" = x; then
      AC_MSG_RESULT(none needed)
   else
      AC_MSG_RESULT($M4FLAGS)
   fi
   M4FFLAGS=`echo $M4FLAGS | sed 's/-s *//g'`
   dnl AC_MSG_NOTICE(M4FFLAGS=$M4FFLAGS)
   AC_SUBST(M4FLAGS)
   AC_SUBST(M4FFLAGS)
])

dnl
dnl Check for an ar(1) utility.
dnl
AC_DEFUN([UD_PROG_AR], [
   dnl AS_MESSAGE([checking for ar utility...])
   case "${AR-unset}" in
       unset) AC_CHECK_PROGS(AR, ar, ar) ;;
       *) AC_CHECK_PROGS(AR, $AR ar, ar) ;;
   esac
   AC_MSG_CHECKING(ar flags)
   case "${ARFLAGS-unset}" in
       unset) ARFLAGS=cru ;;
   esac
   AC_MSG_RESULT($ARFLAGS)
   AC_SUBST(ARFLAGS)
])

dnl
dnl Check for an nm(1) utility.
dnl
AC_DEFUN([UD_PROG_NM], [
   dnl AS_MESSAGE([checking for nm utility...])
   case "${NM-unset}" in
       unset) AC_CHECK_PROGS(NM, nm, nm) ;;
       *) AC_CHECK_PROGS(NM, $NM nm, nm) ;;
   esac
   AC_MSG_CHECKING(nm flags)
   case "${NMFLAGS-unset}" in
       unset) NMFLAGS= ;;
   esac
   AC_MSG_RESULT($NMFLAGS)
   AC_SUBST(NMFLAGS)
])

dnl
dnl Set the top-level source-directory.
dnl
AC_DEFUN([UD_SRCDIR], [
   AC_MSG_CHECKING(for top-level source-directory)
   SRCDIR=`(cd $srcdir && pwd)`
   AC_MSG_RESULT($SRCDIR)
   AC_SUBST(SRCDIR)
])

dnl
dnl UD_CHECK_IEEE
dnl If the 'double' is not an IEEE double
dnl or the 'float' is not and IEEE single,
dnl define NO_IEEE_FLOAT
dnl
AC_DEFUN([UD_CHECK_IEEE], [
   AC_MSG_CHECKING(for IEEE floating point format)
   AC_RUN_IFELSE([AC_LANG_PROGRAM([[
      #ifndef NO_FLOAT_H
      #include <float.h>
      #endif
      #define EXIT_NOTIEEE   1
      #define EXIT_MAYBEIEEE 0]], [[
#if defined(FLT_RADIX) && FLT_RADIX != 2
      return EXIT_NOTIEEE;
#elif defined(DBL_MAX_EXP) && DBL_MAX_EXP != 1024
      return EXIT_NOTIEEE;
#elif defined(DBL_MANT_DIG) && DBL_MANT_DIG != 53
      return EXIT_NOTIEEE;
#elif defined(FLT_MAX_EXP) && !(FLT_MAX_EXP == 1024 || FLT_MAX_EXP == 128)
      return EXIT_NOTIEEE;
#elif defined(FLT_MANT_DIG) && !(FLT_MANT_DIG == 53 || FLT_MANT_DIG == 24)
      return EXIT_NOTIEEE;
#else
      /* (assuming eight bit char) */
      if (sizeof(double) != 8)
          return EXIT_NOTIEEE;
      if (!(sizeof(float) == 4 || sizeof(float) == 8))
          return EXIT_NOTIEEE;

      return EXIT_MAYBEIEEE;
#endif
      ]])],[ac_cv_c_ieeefloat=yes],[ac_cv_c_ieeefloat=no],[:])
   AC_MSG_RESULT($ac_cv_c_ieeefloat)
   if test x$ac_cv_c_ieeefloat = xno; then
      AC_DEFINE(NO_IEEE_FLOAT)
   fi
])

dnl Setup for making a manual-page database.
dnl
AC_DEFUN([UD_MAKEWHATIS],
[
    #
    # NB: We always want to define WHATIS to prevent the
    # $(MANDIR)/$(WHATIS) make(1) target from being just $(MANDIR)/ and
    # conflicting with the (directory creation) target with the same name.
    #
    WHATIS=whatis
    case `uname -sr` in
        BSD/OS*|FreeBSD*)
            # Can't generate a user-database -- only /usr/share/man/whatis.db.
            MAKEWHATIS_CMD=
            ;;
        'IRIX64 6.5'|'IRIX 6.5')
            MAKEWHATIS_CMD='/usr/lib/makewhatis -M $(MANDIR) $(MANDIR)/whatis'
            ;;
        'IRIX 6'*)
            # Can't generate a user-database.
            MAKEWHATIS_CMD=
            ;;
        HP-UX*)
            # Can't generate a user-database -- only /usr/lib/whatis.
            MAKEWHATIS_CMD=
            ;;
        'Linux '*)
            # /usr/sbin/makewhatis doesn't work
            MAKEWHATIS_CMD=
            ;;
        ULTRIX*)
            # Can't generate a user-database -- only /usr/lib/whatis.
            MAKEWHATIS_CMD=
            ;;
        *)
            if test -r /usr/man/windex; then
                WHATIS=windex
            fi
            AC_CHECK_PROGS(prog, catman makewhatis /usr/lib/makewhatis)
            case "$prog" in
                *catman*)
                    MAKEWHATIS_CMD=$prog' -w -M $(MANDIR)'
                    ;;
                *makewhatis*)
                    MAKEWHATIS_CMD=$prog' $(MANDIR)'
                    ;;
            esac
            ;;
    esac
    AC_SUBST(WHATIS)
    AC_SUBST(MAKEWHATIS_CMD)
    AC_MSG_CHECKING(for manual-page index command)
    AC_MSG_RESULT($MAKEWHATIS_CMD)
])


dnl Check for the math library.
dnl
AC_DEFUN([UD_CHECK_LIB_MATH], [
   dnl AS_MESSAGE([checking for math library...])
   case "${MATHLIB}" in
       '')
           AC_CHECK_LIB(c, floor, MATHLIB=,
                   AC_CHECK_LIB(m, floor, MATHLIB=-lm, MATHLIB=))
           ;;
      *)
           AC_MSG_RESULT($MATHLIB (user defined))
           ;;
   esac
   AC_SUBST(MATHLIB)
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
AC_DEFUN([AX_C_FLOAT_WORDS_BIGENDIAN], [
   AC_CACHE_CHECK([whether float word ordering is bigendian],
                  [ax_cv_c_float_words_bigendian], [
      ax_cv_c_float_words_bigendian=unknown
      AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
         double d = 90904234967036810337470478905505011476211692735615632014797120844053488865816695273723469097858056257517020191247487429516932130503560650002327564517570778480236724525140520121371739201496540132640109977779420565776568942592.0;
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
         fi]
      )]
   )
   case $ax_cv_c_float_words_bigendian in
      yes)
         m4_default([$1],[
            AC_DEFINE([FLOAT_WORDS_BIGENDIAN], 1,
                      [Define to 1 if your system stores words within floats
                       with the most significant word first])]) ;;
      no)
         $2 ;;
      *)
         m4_default([$3],
            [AC_MSG_ERROR([Unknown float word ordering. You need to manually preset
                           ax_cv_c_float_words_bigendian=no (or yes) according to your system.])

            ]) ;;
   esac
])# AX_C_FLOAT_WORDS_BIGENDIAN

dnl Find the full path of a header file
dnl
dnl UD_CHECK_HEADER_PATH(file, [action-if-found], [action-if-not-found])
dnl Example:
dnl UD_CHECK_HEADER_PATH([math.h])
dnl AC_MSG_NOTICE([ac_cv_header_path_math_h=$ac_cv_header_path_math_h])
dnl
dnl
AC_DEFUN([UD_CHECK_HEADER_PATH], [
   AS_VAR_PUSHDEF([ac_Path], [ac_cv_header_path_$1])dnl
   AC_CACHE_CHECK(
      [for full path of header file \"$1\"], [ac_Path],
      [AC_PREPROC_IFELSE(
          [AC_LANG_PROGRAM([[#include <$1>]])],
          [AS_VAR_SET([ac_Path], [`sed -n '/\.h"/s/.*"\(.*\)".*/\1/p' conftest.i | grep -m 1 $1`])],
          [AC_MSG_RESULT([not found])]
      )])
   AS_VAR_SET_IF([ac_Path], [$2], [$3])
   dnl eval AS_TR_SH([ac_cv_header_path_$1])=$ac_Path
   AS_VAR_POPDEF([ac_Path])dnl
])

dnl
dnl Check how sed command handling in-place option -i and define SED_I
dnl
AC_DEFUN([UD_PROG_SED_I], [
   AC_REQUIRE([AC_PROG_SED])
   AC_CACHE_CHECK([for sed handling option -i ], ac_cv_SED_I,[
   cat > conftest.sed_i <<EOF
   test str1
EOF
   ac_cv_err=`$SED -i '' -e 's|str1|str2|g' conftest.sed_i 2>&1`
   if test "x$ac_cv_err" = x ; then
      ac_cv_SED_I="$SED -i ''"
   else
      ac_cv_err=`sed -i'' -e 's|str1|str2|g' conftest.sed_i 2>&1`
      if test "x$ac_cv_err" = x ; then
         ac_cv_SED_I="$SED -i''"
      else
         AC_MSG_ERROR("No proper sed -i option found")
      fi
   fi
   AS_UNSET(ac_cv_err)])
   SED_I="$ac_cv_SED_I"
   AC_SUBST(SED_I)
   ${RM} -f conftest.sed_i
])

# LT_AC_CHECK_SHLIB(LIBRARY, FUNCTION,
#                   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#                   [OTHER-LIBRARIES])
# -----------------------------------------------------------
#
# Use a cache variable name containing both the library and function name,
# because the test really is for library $1 defining function $2, not
# just for library $1. Separate tests with the same $1 and different $2s
# may have different results.
#
# Note that using directly AS_VAR_PUSHDEF([ac_Lib], [lt_ac_cv_shlib_$1_$2])
# is asking for troubles, since LT_AC_CHECK_SHLIB($lib, fun) would give
# lt_ac_cv_shlib_$lib_fun, which is definitely not what was meant. Hence
# the AS_LITERAL_IF indirection.
#
# FIXME: This macro is extremely suspicious. It DEFINEs unconditionally,
# whatever the FUNCTION, in addition to not being a *S macro.  Note
# that the cache does depend upon the function we are looking for.
#
AC_DEFUN([LT_AC_CHECK_SHLIB],[
   m4_ifval([$3], , [LT_AH_CHECK_SHLIB([$1])])dnl
   AS_LITERAL_IF([$1],
                 [AS_VAR_PUSHDEF([ac_Lib], [lt_ac_cv_shlib_$1_$2])],
                 [AS_VAR_PUSHDEF([ac_Lib], [lt_ac_cv_shlib_$1''_$2])])dnl
   AC_CACHE_CHECK([for $2 in shared version of -l$1], ac_Lib, [
      lt_ac_check_shlib_save_LIBS=$LIBS
      LIBS="-l$1 $5 $LIBS"
      LT_AC_LINK_SHLIB_IFELSE([AC_LANG_CALL([], [$2])],
         [AS_VAR_SET(ac_Lib, yes)], [AS_VAR_SET(ac_Lib, no)])
      LIBS=$lt_ac_check_shlib_save_LIBS
   ])
   AS_IF([test AS_VAR_GET(ac_Lib) = yes],
      [m4_default([$3], [AC_DEFINE_UNQUOTED(AS_TR_CPP(HAVE_SHLIB$1))])],
      [$4])dnl
   AS_VAR_POPDEF([ac_Lib])dnl
])# AC_CHECK_LIB

# LT_AH_CHECK_SHLIB(LIBNAME)
# ---------------------
m4_define([LT_AH_CHECK_SHLIB],[
   AH_TEMPLATE(AS_TR_CPP(HAVE_SHLIB$1),
      [Define to 1 if you have a shared version of the `]$1[' library (-l]$1[).])])


# LT_AC_LINK_SHLIB_IFELSE(LIBRARYCODE, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------------------
# Try to link LIBRARYCODE into a libtool library.
AC_DEFUN([LT_AC_LINK_SHLIB_IFELSE],[
   m4_ifvaln([$1], [AC_LANG_CONFTEST([$1])])dnl
   rm -rf $objdir
   rm -f conftest.$ac_objext conftest.la
   ac_ltcompile='./libtool --mode=compile $CC -c $CFLAGS $CPPFLAGS conftest.$ac_ext -o conftest.lo >&AS_MESSAGE_LOG_FD'
   ac_ltlink_la='./libtool --mode=link $CC -rpath `pwd` $CFLAGS $LDFLAGS -o libconftest.la conftest.lo $LIBS >&AS_MESSAGE_LOG_FD'
   AS_IF([AC_TRY_EVAL([ac_ltcompile]) &&
          AC_TRY_EVAL([ac_ltlink_la]) &&
          AC_TRY_COMMAND([test -s libconftest.la])],
         [$2],
         [_AC_MSG_LOG_CONFTEST
         m4_ifvaln([$3], [$3])])
   rm -rf $objdir
   rm -f conftest* libconftest*[]dnl
])# _AC_LINK_IFELSE

