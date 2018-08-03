#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c .f .F .f90

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS  = -I$(top_srcdir)/src/include
AM_CPPFLAGS += -I$(srcdir)/../common
AM_CPPFLAGS += -I$(top_builddir)/src/include
LDADD = ${top_builddir}/src/libs/libpnetcdf.la ../common/libtestutils.la -lm


if DECL_MPI_OFFSET
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   # AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
   AM_CFLAGS   = -DHAVE_DECL_MPI_OFFSET
endif

NCMPIGEN  = $(top_builddir)/src/utils/ncmpigen/ncmpigen
NCMPIDIFF = $(top_builddir)/src/utils/ncmpidiff/ncmpidiff

TESTPROGRAMS = simple_xy \
               put_get_all_kinds \
               interoperability_rd \
               interoperability_wr \
               noclobber

M4_SRCS  = put_get_all_kinds.m4 interoperability_rd.m4 interoperability_wr.m4

nodist_put_get_all_kinds_SOURCES = put_get_all_kinds.c
nodist_interoperability_rd_SOURCES = interoperability_rd.c
nodist_interoperability_wr_SOURCES = interoperability_wr.c

M4FLAGS += -I${top_srcdir}/m4

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

check_PROGRAMS = $(TESTPROGRAMS) \
                 compressed

# autimake 1.11.3 has not yet implemented AM_TESTS_ENVIRONMENT
# For newer versions, we can use AM_TESTS_ENVIRONMENT instead
# AM_TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
# AM_TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
# AM_TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;
TESTS_ENVIRONMENT  = TESTPROGRAMS="$(TESTPROGRAMS)" ; export TESTPROGRAMS;
TESTS_ENVIRONMENT += TESTSEQRUN="$(TESTSEQRUN)" ; export TESTSEQRUN;
TESTS_ENVIRONMENT += TESTOUTDIR="$(TESTOUTDIR)" ; export TESTOUTDIR;

TESTS = $(TESTPROGRAMS) seq_runs.sh
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/wrap_runs.sh
SH_LOG_COMPILER =

NC_FILES = $(TESTPROGRAMS:%=$(TESTOUTDIR)/%.nc)

CLEANFILES = $(M4_SRCS:.m4=.c) \
             $(TESTOUTDIR)/put_get_all_kinds.nc.cdf4 \
             $(NC_FILES)

EXTRA_DIST = $(M4_SRCS) seq_runs.sh wrap_runs.sh # redef-good.ncdump

../common/libtestutils.la:
	set -e; cd ../common && $(MAKE) $(MFLAGS) tests

# Some of these tests are designed to run on one process,
# Run them on 4 processes to see if they can handle well
# Some of these tests are designed to run on 4 processes,
# Run them on 2, 4, and 6 processes to see if they can handle well

TESTMPIRUN2  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/2/g'`
TESTMPIRUN4  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/4/g'`
TESTMPIRUN6  = `echo $(TESTMPIRUN) | $(SED) -e 's/NP/6/g'`

ptest ptest4: $(TESTPROGRAMS)
	export PNETCDF_SAFE_MODE=0 && $(TESTMPIRUN4) ./compressed $(srcdir)
	export PNETCDF_SAFE_MODE=1 && $(TESTMPIRUN4) ./compressed $(srcdir)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN4) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	) ; done ; } ; done

ptest2: $(TESTPROGRAMS)
	export PNETCDF_SAFE_MODE=0 && $(TESTMPIRUN2) ./compressed $(srcdir)
	export PNETCDF_SAFE_MODE=1 && $(TESTMPIRUN2) ./compressed $(srcdir)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN2) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	) ; done ; } ; done

ptest6: $(TESTPROGRAMS)
	export PNETCDF_SAFE_MODE=0 && $(TESTMPIRUN6) ./compressed $(srcdir)
	export PNETCDF_SAFE_MODE=1 && $(TESTMPIRUN6) ./compressed $(srcdir)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN6) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	) ; done ; } ; done

ptests: ptest2 ptest4 ptest6
ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

