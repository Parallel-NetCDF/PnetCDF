#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

srcdir = @srcdir@
VPATH = @srcdir@

include ../../macros.make

NC_ZIP_FILES = test_double.nc.gz \
               test_float.nc.gz \
               test_int_cdf5.nc.gz \
               test_int.nc.gz

NC_FILES = $(NC_ZIP_FILES:.nc.gz=.nc)

GARBAGE = $(NC_FILES)

PACKING_LIST = $(NC_ZIP_FILES) Makefile.in

all: $(NC_FILES)

test_double.nc: test_double.nc.gz
	gzip -dc $< > $@

test_float.nc: test_float.nc.gz
	gzip -dc $< > $@

test_int_cdf5.nc: test_int_cdf5.nc.gz
	gzip -dc $< > $@

test_int.nc: test_int.nc.gz
	gzip -dc $< > $@

unzip:
	set -e; for i in $(NC_ZIP_FILES); do ( if [ ! -f $$i ] ; then gzip -dc $(srcdir)/$$i.gz > $$i ; fi ) ; done

testing check verbose_testing:

ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10:


include $(srcdir)/../../rules.make

