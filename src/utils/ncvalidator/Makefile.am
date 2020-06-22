#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(top_builddir)/src/include

if PNETCDF_DEBUG
   AM_CPPFLAGS += -DPNETCDF_DEBUG
endif

bin_PROGRAMS = ncvalidator
ncvalidator_SOURCES = ncvalidator.c

dist_man_MANS = ncvalidator.1

check_PROGRAMS = $(bin_PROGRAMS) tst_open

tst_open_LDADD  = $(top_builddir)/src/libs/libpnetcdf.la
tst_open_LDADD += @NETCDF4_LDFLAGS@ @ADIOS_LDFLAGS@ @NETCDF4_LIBS@ @ADIOS_LIBS@

$(top_builddir)/src/libs/libpnetcdf.la:
	set -e; cd $(top_builddir)/src/libs && $(MAKE) $(MFLAGS)

ENULLPAD_FILES = bad_padding_name_dim.nc1 \
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
                 bad_padding_value_vattr.nc5

EMAXVARS_FILES = bad_nvars.nc1 \
                 bad_nvars.nc2 \
                 bad_nvars.nc5

EUNLIMIT_FILES = bad_unlimited.nc1 \
                 bad_unlimited.nc2 \
                 bad_unlimited.nc5

ENOTNC_FILES =   bad_magic.nc \
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

EVARSIZE_FILES = bad_large_fixed_var.nc2 \
                 bad_large_rec_2_vars.nc2 \
                 bad_large_rec_var.nc2

TST_HDF5_FILES = pad_superblock.h5

ncvalidator$(EXEEXT): ncvalidator.c Makefile
	$(SEQ_CC) -o $@ $<

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
# AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
# AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT += ENULLPAD_FILES="$(ENULLPAD_FILES)" ; export ENULLPAD_FILES;
TESTS_ENVIRONMENT += EMAXVARS_FILES="$(EMAXVARS_FILES)" ; export EMAXVARS_FILES;
TESTS_ENVIRONMENT += EUNLIMIT_FILES="$(EUNLIMIT_FILES)" ; export EUNLIMIT_FILES;
TESTS_ENVIRONMENT += ENOTNC_FILES="$(ENOTNC_FILES)" ; export ENOTNC_FILES;
TESTS_ENVIRONMENT += EVARSIZE_FILES="$(EVARSIZE_FILES)" ; export EVARSIZE_FILES;
TESTS_ENVIRONMENT += TST_HDF5_FILES="$(TST_HDF5_FILES)" ; export TST_HDF5_FILES;

TESTS = xfail_runs.sh wrap_runs.sh
XFAIL_TESTS = xfail_runs.sh

EXTRA_DIST = xfail_runs.sh wrap_runs.sh \
             $(ENULLPAD_FILES) $(EMAXVARS_FILES) \
             $(EUNLIMIT_FILES) $(ENOTNC_FILES) $(EVARSIZE_FILES) \
             $(TST_HDF5_FILES)

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

dist-hook:
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g"            $(distdir)/ncvalidator.c
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/ncvalidator.c
	$(SED_I) -e "s|PNETCDF_RELEASE_VERSION|$(PNETCDF_VERSION)|g"            $(distdir)/ncvalidator.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE_FULL|@PNETCDF_RELEASE_DATE_FULL@|g" $(distdir)/ncvalidator.1
	$(SED_I) -e "s|PNETCDF_RELEASE_DATE|@PNETCDF_RELEASE_DATE@|g"           $(distdir)/ncvalidator.1

