#
# Copyright (C) 2015, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

SUFFIXES = .o .c

CC = $(SEQ_CC)

bin_PROGRAMS = ncoffsets
ncoffsets_SOURCES = ncoffsets.c

dist_man_MANS = ncoffsets.1

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

