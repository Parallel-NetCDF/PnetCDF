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
AM_FFLAGS    = -I$(top_builddir)/src/binding/f77
AM_FCFLAGS   = $(FC_MODINC)$(top_builddir)/src/binding/f90
AM_FCFLAGS  += $(FC_MODINC)$(srcdir)/../common
AM_FCFLAGS  += -I$(top_builddir)/src/binding/f77
LDADD = ${top_builddir}/src/libs/libpnetcdf.la ../common/libtestutils.la -lm

if SIZEOF_MPI_AINT_IS_4
   AM_FFLAGS += $(FC_DEFINE)SIZEOF_MPI_AINT_IS_4
endif

if DECL_MPI_OFFSET
   # Do not add to AM_CPPFLAGS, as it will also be used by Fortran programs
   # AM_CPPFLAGS += -DHAVE_DECL_MPI_OFFSET
   AM_CFLAGS   = -DHAVE_DECL_MPI_OFFSET
   AM_FFLAGS  += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
   AM_FCFLAGS += $(FC_DEFINE)HAVE_DECL_MPI_OFFSET
endif

NCMPIGEN  = $(top_builddir)/src/utils/ncmpigen/ncmpigen
NCMPIDIFF = $(top_builddir)/src/utils/ncmpidiff/ncmpidiff

TESTPROGRAMS = ncmpi_vars_null_stride \
               vectors \
               collective_error \
               test_varm \
               alignment_test \
               flexible \
               flexible2 \
               flexible_varm \
               nonblocking \
               noclobber \
               record \
               inq_num_vars \
               varn_int \
               modes \
               one_record \
               inq_recsize \
               test_vard \
               test_vard_rec \
               test_vard_multiple \
               varn_contig \
               ivarn \
               check_striping \
               add_var \
               buftype_free \
               last_large_var \
               check_type \
               test_erange \
               scalar \
               null_args \
               tst_dimsizes \
               mix_collectives \
               large_var_cdf5 \
               tst_max_var_dims \
               tst_info \
               tst_vars_fill \
               tst_def_var_fill

M4_SRCS  = put_all_kinds.m4 erange_fill.m4 tst_vars_fill.m4 iput_all_kinds.m4

nodist_put_all_kinds_SOURCES = put_all_kinds.c
nodist_erange_fill_SOURCES = erange_fill.c
nodist_tst_vars_fill_SOURCES = tst_vars_fill.c
nodist_iput_all_kinds_SOURCES = iput_all_kinds.c

if ENABLE_ERANGE_FILL
   M4FLAGS += -DERANGE_FILL
   TESTPROGRAMS += erange_fill
endif

M4FLAGS += -I${top_srcdir}/m4

$(M4_SRCS:.m4=.c): Makefile

.m4.c:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

if HAS_FORTRAN
   TESTPROGRAMS += varn_intf \
                   attrf \
                   buftype_freef \
                   put_parameter \
                   test_vardf

       varn_intf_SOURCES = varn_intf.f
           attrf_SOURCES = attrf.f
   buftype_freef_SOURCES = buftype_freef.f
   put_parameter_SOURCES = put_parameter.f
      test_vardf_SOURCES = test_vardf.F
if HAVE_MPI_MOD
   TESTPROGRAMS += inq_num_varsf \
                   inq_recsizef \
                   test_vardf90 \
                   varn_real

   inq_num_varsf_SOURCES = inq_num_varsf.f90
    inq_recsizef_SOURCES = inq_recsizef.f90
    test_vardf90_SOURCES = test_vardf90.f90
       varn_real_SOURCES = varn_real.f90
endif
endif

check_PROGRAMS = $(TESTPROGRAMS) put_all_kinds redef1 tst_vars_fill \
                 iput_all_kinds

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

if BUILD_DRIVER_DW
TESTDW = 1
TESTS_ENVIRONMENT += TESTDW=1 ; export TESTDW;
NC_FILES += $(TESTOUTDIR)/put_all_kinds_dw.nc.cdf1 \
            $(TESTOUTDIR)/put_all_kinds_dw.nc.cdf2 \
            $(TESTOUTDIR)/put_all_kinds_dw.nc.cdf5 \
            $(TESTOUTDIR)/iput_all_kinds_dw.nc.cdf1 \
            $(TESTOUTDIR)/iput_all_kinds_dw.nc.cdf2 \
            $(TESTOUTDIR)/iput_all_kinds_dw.nc.cdf5
endif

BURST_BUFFER_FILES = $(NC_FILES:=_*.meta) $(NC_FILES:=_*.data)

CLEANFILES = $(M4_SRCS:.m4=.c) core core.* *.gcda *.gcno *.gcov gmon.out \
             $(TESTOUTDIR)/testfile.nc \
             $(TESTOUTDIR)/redef1.nc \
             $(TESTOUTDIR)/put_all_kinds.nc.cdf1 \
             $(TESTOUTDIR)/put_all_kinds.nc.cdf2 \
             $(TESTOUTDIR)/put_all_kinds.nc.cdf5 \
             $(TESTOUTDIR)/iput_all_kinds.nc.cdf1 \
             $(TESTOUTDIR)/iput_all_kinds.nc.cdf2 \
             $(TESTOUTDIR)/iput_all_kinds.nc.cdf5 \
             $(NC_FILES)

EXTRA_DIST = $(M4_SRCS) seq_runs.sh redef-good.ncdump \
             wrap_runs.sh

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
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN4) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	if [ -n "$(TESTDW)" ]; then ( \
	export PNETCDF_HINTS="nc_dw=enable;nc_dw_dirname=${TESTOUTDIR};nc_dw_overwrite=enable" ; \
	$(TESTMPIRUN4) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	unset PNETCDF_HINTS ; \
	) ; fi ; \
	) ; done ; } ; done

ptest2: $(TESTPROGRAMS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN2) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	if [ -n "$(TESTDW)" ]; then ( \
	export PNETCDF_HINTS="nc_dw=enable;nc_dw_dirname=${TESTOUTDIR};nc_dw_overwrite=enable" ; \
	$(TESTMPIRUN2) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	unset PNETCDF_HINTS ; \
	) ; fi ; \
	) ; done ; } ; done

ptest6: $(TESTPROGRAMS)
	for j in 0 1 ; do { \
	export PNETCDF_SAFE_MODE=$$j ; \
	set -e; for i in $(TESTPROGRAMS); do ( \
	$(TESTMPIRUN_6) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	if [ -n "$(TESTDW)" ]; then ( \
	export PNETCDF_HINTS="nc_dw=enable;nc_dw_dirname=${TESTOUTDIR};nc_dw_overwrite=enable" ; \
	$(TESTMPIRUN6) ./$$i $(TESTOUTDIR)/$$i.nc ; \
	unset PNETCDF_HINTS ; \
	) ; fi ; \
	) ; done ; } ; done

ptests: ptest2 ptest4 ptest6
ptest8 ptest10:

# build check targets but not invoke
tests-local: all $(check_PROGRAMS)

.PHONY: ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10

