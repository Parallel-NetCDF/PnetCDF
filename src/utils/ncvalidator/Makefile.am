#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

if PNETCDF_DEBUG
   AM_CPPFLAGS = -DPNETCDF_DEBUG
endif

bin_PROGRAMS = ncvalidator
ncvalidator_SOURCES = ncvalidator.c

dist_man_MANS = ncvalidator.1

check_PROGRAMS = $(bin_PROGRAMS)

BAD_FILES = bad_magic.nc \
            bad_nvars.nc1 \
            bad_nvars.nc2 \
            bad_nvars.nc5 \
            bad_padding_name_dim.nc1 \
            bad_padding_name_dim.nc2 \
            bad_padding_name_dim.nc5 \
            bad_padding_name_gattr.nc1 \
            bad_padding_name_gattr.nc2 \
            bad_padding_name_gattr.nc5 \
            bad_padding_name_var.nc1 \
            bad_padding_name_var.nc2 \
            bad_padding_name_var.nc5 \
            bad_padding_name_vattr.nc1 \
            bad_padding_name_vattr.nc2 \
            bad_padding_name_vattr.nc5 \
            bad_padding_value_gattr.nc1 \
            bad_padding_value_gattr.nc2 \
            bad_padding_value_gattr.nc5 \
            bad_padding_value_vattr.nc1 \
            bad_padding_value_vattr.nc2 \
            bad_padding_value_vattr.nc5 \
            bad_tag_dimension.nc1 \
            bad_tag_dimension.nc2 \
            bad_tag_dimension.nc5 \
            bad_tag_gattr.nc1 \
            bad_tag_gattr.nc2 \
            bad_tag_gattr.nc5 \
            bad_tag_vattr.nc1 \
            bad_tag_vattr.nc2 \
            bad_tag_vattr.nc5 \
            bad_tag_variable.nc1 \
            bad_tag_variable.nc2 \
            bad_tag_variable.nc5

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
# AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
# AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT += BAD_FILES="$(BAD_FILES)" ; export BAD_FILES;

TESTS = xfail_runs.sh
XFAIL_TESTS = xfail_runs.sh

EXTRA_DIST = xfail_runs.sh $(BAD_FILES)

CLEANFILES = *.gcno

tests-local: all

dist-hook:
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g" $(distdir)/ncvalidator.c
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g" $(distdir)/ncvalidator.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE2|`date '+%Y-%m-%d'`|g"   $(distdir)/ncvalidator.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|`date '+%e %b %Y'`|g"    $(distdir)/ncvalidator.1

